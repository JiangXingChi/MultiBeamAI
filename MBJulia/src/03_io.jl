# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  I/O — 数据加载与结果输出                                                     ║
# ║                                                                              ║
# ║  LoadBathymetry: 从 XYZ 文本文件读取 lon/lat/depth 三列                       ║
# ║  WriteGridOutput: 将过滤后的格网化结果写入 XYZ 文件                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ── 加载原始波束点云 ──
#
# 目的：从空格分隔的 XYZ 文件读取全部波束点，拆为三个独立向量。
# 为什么拆成三列而非保留矩阵：后续各阶段只用部分列（如 depths），
# 拆开避免重复切片开销，且语义更清晰。
function LoadBathymetry(filepath::String)
    println("[1/6] 加载原始波束数据...")
    @time raw = readdlm(filepath, Float64; comments=false)
    lons   = raw[:, 1]
    lats   = raw[:, 2]
    depths = raw[:, 3]
    n = length(depths)
    println("  已加载 $n 个波束点")
    return lons, lats, depths, n
end


# ── 输出格网化 XYZ ──
#
# 目的：遍历所有非空格网，取保留波束点的均值深度和格网中心坐标，写入文件。
# 输出格式：lon_center lat_center mean_depth（空格分隔，标准 XYZ）
#
# 格网中心坐标公式：第 col 列中心 = lon_min + (col - 0.5) × lon_step
# 因为第 1 列覆盖 [lon_min, lon_min + lon_step)，中心在半步长处。
function WriteGridOutput(filepath::String, cell_to_indices, keep_mask,
                          depths, geo::GridGeometry)
    println("\n[6/6] 输出格网化结果...")

    cells_output = 0
    cells_empty  = 0

    @time open(filepath, "w") do io
        for key in keys(cell_to_indices)
            idxs = cell_to_indices[key]
            kept_idxs = idxs[keep_mask[idxs]]
            if isempty(kept_idxs)
                cells_empty += 1
                continue
            end

            cells_output += 1
            col, row = key
            lon_c = geo.lon_min + (col - 0.5) * geo.lon_step
            lat_c = geo.lat_min + (row - 0.5) * geo.lat_step
            depth_mean = mean(depths[kept_idxs])
            @printf(io, "%.8f %.8f %.3f\n", lon_c, lat_c, depth_mean)
        end
    end

    out_size = stat(filepath).size / 1e6
    @printf("  已输出: %s（%.1f MB）\n", filepath, out_size)
    @printf("  输出格网: %d，全噪声被清空格网: %d\n", cells_output, cells_empty)
    return cells_output, cells_empty
end
