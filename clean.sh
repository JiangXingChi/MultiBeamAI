#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════════════╗
# ║  多波束处理管线 — 清理脚本                                            ║
# ║                                                                      ║
# ║  用法：                                                               ║
# ║    bash clean.sh all      清理全部 (MBSystem + Julia结果 + R结果)    ║
# ║    bash clean.sh mb       删除整个 MBSystem/ (s7k/esf/result)       ║
# ║    bash clean.sh julia    清空 MBJulia/result/                      ║
# ║    bash clean.sh r        清空 MBR/result/                          ║
# ║    bash clean.sh output   仅清理最终输出文件 (保留中间产物)          ║
# ║    bash clean.sh result   清空 Result/ (归档输出)                    ║
# ╚══════════════════════════════════════════════════════════════════════╝

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
MB_SYSTEM="$PROJECT_ROOT/MBSystem"
JULIA_RESULT="$PROJECT_ROOT/MBJulia/result"
R_RESULT="$PROJECT_ROOT/MBR/result"
RESULT_DIR="$PROJECT_ROOT/Result"

clear_dir() {
    local dir="$1"; local label="$2"
    if [ -d "$dir" ]; then
        rm -rf "$dir"/*
        echo "  ✓ $label 已清空"
    else
        mkdir -p "$dir"
        echo "  ✓ $label 已创建（原先不存在）"
    fi
}

rebuild_mb() {
    echo "  → 删除 $MB_SYSTEM"
    rm -rf "$MB_SYSTEM"
    mkdir -p "$MB_SYSTEM/data" "$MB_SYSTEM/result"
    echo "  ✓ MBSystem/ 已重建 (data/ + result/)"
}

clean_output() {
    echo "  → 清理最终输出 (3 txt + 3 png)"
    rm -f "$JULIA_RESULT"/gridded_1m_xyz.txt
    rm -f "$JULIA_RESULT"/gridded_1m_clean.txt
    rm -f "$JULIA_RESULT"/gridded_1m_interpolated.txt
    rm -f "$JULIA_RESULT"/boundary_1m.txt
    rm -f "$R_RESULT"/gridded_1m_xyz.png
    rm -f "$R_RESULT"/gridded_1m_clean.png
    rm -f "$R_RESULT"/gridded_1m_interpolated.png
    echo "  ✓ 最终输出已清理"
}

clean_result() {
    if [ -d "$RESULT_DIR" ]; then
        rm -rf "$RESULT_DIR"/*
        echo "  ✓ Result/ 已清空"
    else
        echo "  Result/ 不存在，跳过"
    fi
}

MODE="${1:-all}"

echo "========================================"
echo "多波束管线清理"
echo "========================================"

case "$MODE" in
    all)
        rebuild_mb
        clear_dir "$JULIA_RESULT" "MBJulia/result/"
        clear_dir "$R_RESULT"  "MBR/result/"
        clean_result
        ;;
    mb)
        rebuild_mb
        ;;
    julia)
        clear_dir "$JULIA_RESULT" "MBJulia/result/"
        ;;
    r)
        clear_dir "$R_RESULT"  "MBR/result/"
        ;;
    output)
        clean_output
        ;;
    result)
        clean_result
        ;;
    *)
        echo "未知模式: $MODE"
        echo "可用: all | mb | julia | r | output | result"
        exit 1
        ;;
esac

echo ""
echo "完成"
