# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Pipeline — 编排管线                                                          ║
# ║                                                                              ║
# ║  Run 是唯一的公共入口，按顺序调用六个处理阶段。                                  ║
# ║  输入输出路径为必传参数，由 run.jl 统一配置。                                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

function Run(;
    input_file::String,
    output_file::String,
)
    println("="^60)
    println("两阶段迭代滤波 — MBJulia")
    println("="^60)

    # Stage 1: 加载
    lons, lats, depths, n_total = LoadBathymetry(input_file)

    # Stage 2: 计算格网几何
    geo = ComputeGridGeometry(lons, lats, GRID_RESOLUTION_M)

    # Stage 3: 点→格分箱
    cell_to_indices = BinPointsToCells(lons, lats, geo)

    # Stage 4: 第一阶段 — 格网内深度间隙聚类
    cell_depth, cell_simple, cell_questionable, questionable_keys =
        Stage1GapFilter(cell_to_indices, depths)

    # Stage 5: 第二阶段 — 迭代空间一致性校验
    keep_mask, _ = Stage2SpatialFilter(
        questionable_keys, cell_depth, cell_simple, cell_questionable, n_total)

    # Stage 6: 输出
    cells_output, cells_empty = WriteGridOutput(
        output_file, cell_to_indices, keep_mask, depths, geo)

    # ── 最终汇总 ──
    n_kept = count(keep_mask)
    println("\n" * "="^60)
    println("滤波完成 — 最终汇总")
    println("="^60)
    @printf("  测区尺寸: %.1f m × %.1f m，共 %d 个格网（%d 非空）\n",
            geo.width_m, geo.height_m, geo.cols * geo.rows, length(cell_to_indices))
    @printf("  第一阶段: %d 简单 + %d 含间隙\n",
            length(cell_simple), length(questionable_keys))
    @printf("  第二阶段: 保留 %d / %d（%.1f%%），剔除噪点 %d（%.1f%%）\n",
            n_kept, n_total, n_kept / n_total * 100,
            n_total - n_kept, (n_total - n_kept) / n_total * 100)
    @printf("  格网输出: %d 个格网 → %s\n", cells_output, output_file)
end
