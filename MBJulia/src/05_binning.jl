# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Stage 3 — 点→格分箱                                                        ║
# ║                                                                              ║
# ║  目的：每个波束点按经纬度落入一个 1m×1m 格子，构建 cell_key → [global_indices] 映射。║
# ║  为什么用 Dict 不是二维数组：网格大多为空（湖泊不规则），Dict 只存非空单元，         ║
# ║  节省内存且遍历更快（只遍历有数据的格）。                                        ║
# ║                                                                              ║
# ║  索引从 1 开始（Julia 惯例）。                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

function BinPointsToCells(lons, lats, geo::GridGeometry)
    println("[3/6] 将波束点分配到格网...")

    @time begin
        # (lon - lon_min) / lon_step = 该点在第几列（浮点），向下取整+1=列号
        ix = floor.(Int, (lons .- geo.lon_min) ./ geo.lon_step) .+ 1
        iy = floor.(Int, (lats .- geo.lat_min) ./ geo.lat_step) .+ 1
        # clamp 防止浮点舍入导致的越界（如算出 0 或 cols+1）
        ix = clamp.(ix, 1, geo.cols)
        iy = clamp.(iy, 1, geo.rows)
    end

    # 构建字典：一次遍历，O(n)
    cell_to_indices = Dict{Tuple{Int,Int}, Vector{Int}}()
    for i in 1:length(lons)
        key = (ix[i], iy[i])
        if haskey(cell_to_indices, key)
            push!(cell_to_indices[key], i)
        else
            cell_to_indices[key] = [i]
        end
    end

    n_cells = length(cell_to_indices)
    @printf("  有效（非空）格网数: %d\n", n_cells)
    return cell_to_indices
end
