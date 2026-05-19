#!/usr/bin/env julia
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  S7K 多波束 Julia 处理 — 完整管线                                             ║
# ║                                                                              ║
# ║  阶段 0：简单格网平均  →  gridded_1m_xyz.txt           (文件①)                ║
# ║  阶段 A：深度间隙聚类 + 迭代空间投票  →  gridded_1m_voted.txt                  ║
# ║  阶段 B：P99 全局离群剔除  →  gridded_1m_clean.txt     (文件②)                ║
# ║  阶段 C：湖泊边界检测  →  boundary_1m.txt                                    ║
# ║  阶段 D：迭代 IDW 插值  →  gridded_1m_interpolated.txt (文件③)                ║
# ║                                                                              ║
# ║  优化：XYZ 只加载一次，几何只计算一次，阶段0向量化。                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using Pkg
Pkg.activate(@__DIR__)
using MBJulia
using DelimitedFiles
using Printf
using Statistics
using ProgressMeter

# ═══════════════════════════════════════════════════════════════════════════════
# 配置 — 只改这里
# ═══════════════════════════════════════════════════════════════════════════════

INPUT_FILE = joinpath(@__DIR__, "..", "MBSystem", "result", "bathymetry_xyz.txt")
RESULT_DIR = joinpath(@__DIR__, "result")

OUTPUT_GRID_RAW  = joinpath(RESULT_DIR, "gridded_1m_xyz.txt")           # 文件①
OUTPUT_VOTED     = joinpath(RESULT_DIR, "gridded_1m_voted.txt")
OUTPUT_CLEAN     = joinpath(RESULT_DIR, "gridded_1m_clean.txt")         # 文件②
OUTPUT_ANOMALIES = joinpath(RESULT_DIR, "gridded_1m_anomalies.txt")
OUTPUT_BOUNDARY  = joinpath(RESULT_DIR, "boundary_1m.txt")
OUTPUT_INTERP    = joinpath(RESULT_DIR, "gridded_1m_interpolated.txt")  # 文件③

CUTOFF_QUANTILE = 0.99  # P99 阈值

# ═══════════════════════════════════════════════════════════════════════════════
# 前置：加载 XYZ 一次，计算几何一次（共享给所有阶段）
# ═══════════════════════════════════════════════════════════════════════════════

function load_xyz(filepath::String)
    isfile(filepath) || error("输入文件不存在: $filepath")
    raw = readdlm(filepath, Float64)
    lons, lats, depths = raw[:,1], raw[:,2], raw[:,3]
    @printf("  波束点: %d,  深度: %.2f ~ %.2f m\n", length(depths), minimum(depths), maximum(depths))
    return lons, lats, depths
end

# 使用 MBJulia 的格网几何计算（唯一真相源）
function build_grid_geometry(lons, lats)
    return MBJulia.ComputeGridGeometry(lons, lats, MBJulia.GRID_RESOLUTION_M)
end

# ── 启动 ──
println("="^60)
println("加载数据 + 计算格网几何（全阶段共享）")
println("="^60)

beam_lons, beam_lats, beam_depths = load_xyz(INPUT_FILE)
geo = build_grid_geometry(beam_lons, beam_lats)

# 为后续阶段预计算 col_center / row_center（闭包，避免重复计算）
col_center(c) = geo.lon_min + (c - 0.5) * geo.lon_step
row_center(r) = geo.lat_min + (r - 0.5) * geo.lat_step

mkpath(RESULT_DIR)

# ═══════════════════════════════════════════════════════════════════════════════
# 阶段 0：简单格网平均（向量化版 — 比逐点循环快 ~50x）
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("阶段 0 — 简单格网平均（不做聚类投票）")
println("="^60)

# 一次性向量化计算所有点的行列索引
@time begin
    ix = clamp.(floor.(Int, (beam_lons .- geo.lon_min) ./ geo.lon_step) .+ 1, 1, geo.cols)
    iy = clamp.(floor.(Int, (beam_lats .- geo.lat_min) ./ geo.lat_step) .+ 1, 1, geo.rows)
end

# 累积求和 — 使用 @inbounds 的快速标量循环（Julia 会良好优化）
sum_depth = zeros(Float64, geo.rows, geo.cols)
count_mat = zeros(Int, geo.rows, geo.cols)
@showprogress desc="  格网平均..." for i in eachindex(beam_depths)
    @inbounds begin
        r, c = iy[i], ix[i]
        sum_depth[r, c] += beam_depths[i]
        count_mat[r, c] += 1
    end
end

# 写入输出 — 使用 printf 批量写入
open(OUTPUT_GRID_RAW, "w") do io
    for r in 1:geo.rows, c in 1:geo.cols
        count_mat[r, c] > 0 || continue
        @printf(io, "%.8f\t%.8f\t%.3f\n",
                col_center(c), row_center(r),
                sum_depth[r, c] / count_mat[r, c])
    end
end

n_avg = count(x -> x > 0, count_mat)
@printf("  格网平均: %d 格网 → %s\n", n_avg, OUTPUT_GRID_RAW)

# ═══════════════════════════════════════════════════════════════════════════════
# 阶段 A：深度间隙聚类 + 多次迭代空间一致性校验
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("阶段 A — 深度间隙聚类 + 多次迭代空间一致性校验 (MBJulia)")
println("="^60)
Run(input_file = INPUT_FILE, output_file = OUTPUT_VOTED)

# ═══════════════════════════════════════════════════════════════════════════════
# 阶段 B：全局分位数离群剔除
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("阶段 B — 全局分位数离群剔除")
println("="^60)

data   = readdlm(OUTPUT_VOTED, Float64)
lons   = data[:,1]; lats = data[:,2]; depths = data[:,3]
n_total = length(depths)

@printf("  点数: %d,  水深: %.2f ~ %.2f m\n", n_total, minimum(depths), maximum(depths))
@printf("  均值: %.2f m  中位: %.2f m  σ: %.2f m\n", mean(depths), median(depths), std(depths))

threshold = quantile(depths, CUTOFF_QUANTILE)
@printf("  P%.0f 阈值: %.2f m\n", CUTOFF_QUANTILE * 100, threshold)

good_mask = depths .<= threshold
n_good, n_bad = count(good_mask), n_total - n_good

@printf("  保留: %d (%.2f%%),  剔除: %d (%.2f%%)\n", n_good, 100*n_good/n_total, n_bad, 100*n_bad/n_total)

open(OUTPUT_CLEAN, "w") do io
    writedlm(io, [lons[good_mask] lats[good_mask] depths[good_mask]])
end
@printf("  → %s (%d 行)  ← 文件②\n", OUTPUT_CLEAN, n_good)

if n_bad > 0
    open(OUTPUT_ANOMALIES, "w") do io
        writedlm(io, [lons[.!good_mask] lats[.!good_mask] depths[.!good_mask]])
    end
    @printf("  异常点 → %s (%d 行)\n", OUTPUT_ANOMALIES, n_bad)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 阶段 C：湖泊边界检测
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("阶段 C — 湖泊边界检测")
println("="^60)
include(joinpath(@__DIR__, "src", "09_boundary.jl"))
DetectBoundary(OUTPUT_CLEAN, OUTPUT_BOUNDARY)

# ═══════════════════════════════════════════════════════════════════════════════
# 阶段 D：迭代 IDW 插值
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("阶段 D — 迭代 IDW 插值（格网坐标多边形判定）")
println("="^60)
include(joinpath(@__DIR__, "src", "10_interpolation.jl"))
InterpolateGrid(OUTPUT_CLEAN, OUTPUT_BOUNDARY, OUTPUT_INTERP;
    power=2, min_neighbors=3, initial_radius=5, max_radius=50)

# ═══════════════════════════════════════════════════════════════════════════════
# 清理中间产物
# ═══════════════════════════════════════════════════════════════════════════════

rm(OUTPUT_VOTED; force=true)
rm(OUTPUT_ANOMALIES; force=true)

# ═══════════════════════════════════════════════════════════════════════════════
# 汇总
# ═══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^60)
println("Julia 处理完成 — 三文件汇总")
println("="^60)
println("  文件①  gridded_1m_xyz.txt           $n_avg 格网 — 纯MBSystem+格网平均")
println("  文件②  gridded_1m_clean.txt         $n_good 格网 — 聚类投票+P99")
println("  文件③  gridded_1m_interpolated.txt  轮廓内100%填充 — 聚类投票+P99+IDW插值")
println("="^60)
