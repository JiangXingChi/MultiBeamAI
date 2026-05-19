#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  S7K 多波束批量处理                                                          ║
# ║                                                                              ║
# ║  自动发现 S7K/ 下所有湖泊，逐个处理：                                          ║
# ║    s7k → mbclean → XYZ → Julia(5阶段) → R可视化 → 存档 → 清理 → 下一个         ║
# ║                                                                              ║
# ║  输出：Result/湖泊名/gridded_1m_{xyz,clean,interpolated}.{txt,png}             ║
# ║                                                                              ║
# ║  用法：bash run.sh                                                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

set -euo pipefail

# ═══════════════════════════════════════════════════════════════════════════════
# 配置
# ═══════════════════════════════════════════════════════════════════════════════
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
S7K_BASE="$PROJECT_ROOT/S7K"
MB_DIR="$PROJECT_ROOT/MBSystem"
JULIA_DIR="$PROJECT_ROOT/MBJulia"
R_DIR="$PROJECT_ROOT/MBR"
RESULT_DIR="$PROJECT_ROOT/Result"

MBCLEAN_MODE="1"
MBCLEAN_SIGMA="2.5"
MBCLEAN_SPIKE="3.0/1/1"
MB_FORMAT=88

# 进度条宽度
BAR_WIDTH=30

# ═══════════════════════════════════════════════════════════════════════════════
# 前置检查
# ═══════════════════════════════════════════════════════════════════════════════
command -v mbinfo  >/dev/null 2>&1 || { echo "错误: mbinfo 未找到"; exit 1; }
command -v mbclean >/dev/null 2>&1 || { echo "错误: mbclean 未找到"; exit 1; }
command -v mblist  >/dev/null 2>&1 || { echo "错误: mblist 未找到"; exit 1; }
command -v julia   >/dev/null 2>&1 || { echo "错误: julia 未找到"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "错误: Rscript 未找到"; exit 1; }

if [ ! -d "$S7K_BASE" ]; then
    echo "错误: S7K 目录不存在: $S7K_BASE"; exit 1
fi

# 发现湖泊（S7K 下的子目录）
mapfile -t LAKES < <(find "$S7K_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
if [ ${#LAKES[@]} -eq 0 ]; then
    echo "错误: S7K/ 下没有湖泊子目录"; exit 1
fi

echo "发现 ${#LAKES[@]} 个湖泊: ${LAKES[*]}"
echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# 进度条
# ═══════════════════════════════════════════════════════════════════════════════
draw_progress() {
    local current=$1 total=$2 label=$3 status=${4:-"处理中"}
    local filled=$(( current * BAR_WIDTH / total ))
    local empty=$(( BAR_WIDTH - filled ))
    local bar
    bar=$(printf '█%.0s' $(seq 1 $filled))
    bar+=$(printf '░%.0s' $(seq 1 $empty))
    printf "\r[%s] %d/%d  %-8s  %s" "$bar" "$current" "$total" "$label" "$status"
}

# ═══════════════════════════════════════════════════════════════════════════════
# 单个湖泊处理
# ═══════════════════════════════════════════════════════════════════════════════
process_lake() {
    local lake_name="$1"
    local s7k_dir="$S7K_BASE/$lake_name"
    local s7k_count
    s7k_count=$(find "$s7k_dir" -maxdepth 1 -name '*.s7k' | wc -l)

    if [ "$s7k_count" -eq 0 ]; then
        echo "  ⚠ $lake_name: 无 s7k 文件，跳过"
        return 1
    fi

    # ── Step 1: 准备 MBSystem（直接重建，不用 cd，全部用绝对路径）──
    rm -rf "$MB_DIR"
    mkdir -p "$MB_DIR/data" "$MB_DIR/result"

    # ── Step 2: 复制 s7k ──
    cp "$s7k_dir"/*.s7k "$MB_DIR/data/"

    # ── Step 3: datalist ──
    for f in "$MB_DIR/data"/*.s7k; do echo "$f $MB_FORMAT"; done > "$MB_DIR/datalist.mb-1"

    # ── Step 4: 数据概览 ──
    {
        echo "湖泊: $lake_name"
        echo "s7k 数: $s7k_count"
        echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        mbinfo -I "$MB_DIR/datalist.mb-1" 2>&1
    } > "$MB_DIR/result/aggregate_info.txt"

    # ── Step 5: mbclean ──
    echo "  滤波中..."
    rm -f "$MB_DIR/data"/*.esf "$MB_DIR/data"/*.par
    mbclean -I "$MB_DIR/datalist.mb-1" \
        -M"$MBCLEAN_MODE" -N"$MBCLEAN_SIGMA" -S"$MBCLEAN_SPIKE" -V 2>&1 | grep -E '(flagged|cleaned|records)' || true

    # ── Step 6: 导出 XYZ ──
    echo "  导出 XYZ 点云..."
    mblist -I "$MB_DIR/datalist.mb-1" -FXYZ -MA 2>/dev/null \
      | awk 'BEGIN{c=0} { if ($7 < 0) {printf "%.8f %.8f %.3f\n", $2, $3, -$7; c++} } END{printf "  导出 %d 波束点\n", c >> "/dev/stderr"}' \
      > "$MB_DIR/result/bathymetry_xyz.txt"

    # ── Step 7: Julia ──
    mkdir -p "$JULIA_DIR/result"
    julia --project="$JULIA_DIR" "$JULIA_DIR/run.jl" 2>&1 | grep -E '(阶段|完成|文件|格网平均|保留|插值填充|轮廓内总计|边界|迭代)' || true

    # ── Step 8: R 可视化 ──
    mkdir -p "$R_DIR/result"
    Rscript "$R_DIR/marmap_viz.R" > /dev/null 2>&1
    echo "  可视化完成"

    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# 主循环
# ═══════════════════════════════════════════════════════════════════════════════

START_TIME=$(date +%s)
TOTAL=${#LAKES[@]}

for i in "${!LAKES[@]}"; do
    lake="${LAKES[$i]}"
    current=$((i + 1))

    draw_progress "$current" "$TOTAL" "$lake" "处理中..."

    if process_lake "$lake"; then
        # ── 存档 ──
        mkdir -p "$RESULT_DIR/$lake"
        cp "$JULIA_DIR/result"/gridded_1m_xyz.txt           "$RESULT_DIR/$lake/" 2>/dev/null || true
        cp "$JULIA_DIR/result"/gridded_1m_clean.txt         "$RESULT_DIR/$lake/" 2>/dev/null || true
        cp "$JULIA_DIR/result"/gridded_1m_interpolated.txt  "$RESULT_DIR/$lake/" 2>/dev/null || true
        cp "$R_DIR/result"/gridded_1m_xyz.png               "$RESULT_DIR/$lake/" 2>/dev/null || true
        cp "$R_DIR/result"/gridded_1m_clean.png             "$RESULT_DIR/$lake/" 2>/dev/null || true
        cp "$R_DIR/result"/gridded_1m_interpolated.png      "$RESULT_DIR/$lake/" 2>/dev/null || true

        # ── 清理工作目录（不动 Result/）──
        rm -rf "$MB_DIR"
        mkdir -p "$MB_DIR/data" "$MB_DIR/result"
        find "$JULIA_DIR/result" -mindepth 1 -delete 2>/dev/null || true
        find "$R_DIR/result" -mindepth 1 -delete 2>/dev/null || true

        draw_progress "$current" "$TOTAL" "$lake" "✓ 完成"
        echo ""
    else
        draw_progress "$current" "$TOTAL" "$lake" "✗ 跳过"
        echo ""
    fi
done

# ═══════════════════════════════════════════════════════════════════════════════
ELAPSED=$(( $(date +%s) - START_TIME ))
echo ""
echo "============================================"
echo "全部完成 (耗时 ${ELAPSED}s)"
echo "============================================"
echo ""

# 展示结果
for lake in "${LAKES[@]}"; do
    if [ -d "$RESULT_DIR/$lake" ]; then
        echo "  $lake/"
        for f in "$RESULT_DIR/$lake"/*; do
            [ -f "$f" ] && printf "    %-6s %s\n" "$(du -h "$f" | cut -f1)" "$(basename "$f")"
        done
    fi
done

echo ""
echo "按 Enter 退出..."
read -r
