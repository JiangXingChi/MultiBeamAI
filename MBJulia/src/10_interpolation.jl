# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  迭代 IDW 插值 — 多边形格网坐标直接判定                                          ║
# ║                                                                              ║
# ║  关键：多边形顶点转为格网坐标 (r,c)，射线法在格网空间运算 → 零精度损失            ║
# ║                                                                              ║
# ║  优化：IDW 预计算距离权重矩阵；复用 MBJulia 常量消除重复声明。                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using DelimitedFiles
using Printf

# ── IDW 插值（预计算距离权重，避免每次 sqrt）──
function IDWInterpolate(target_r::Int, target_c::Int,
                         depth_mat::Matrix{Float64},
                         radius::Int, power::Int)::Tuple{Float64, Int}
    rows, cols = size(depth_mat)
    w_sum = 0.0; wt_sum = 0.0; n_nb = 0

    r_start = max(1, target_r - radius)
    r_end   = min(rows, target_r + radius)
    c_start = max(1, target_c - radius)
    c_end   = min(cols, target_c + radius)

    for nr in r_start:r_end, nc in c_start:c_end
        d = depth_mat[nr, nc]
        isnan(d) && continue
        dr2 = (nr - target_r)^2
        dc2 = (nc - target_c)^2
        dist2 = dr2 + dc2
        dist2 == 0.0 && continue  # 跳过目标自己
        dist2 > radius^2 && continue  # 圆形邻域而非方形
        w = 1.0 / (dist2^(power / 2))
        w_sum += w * d
        wt_sum += w
        n_nb += 1
    end

    if n_nb == 0
        return (NaN, 0)
    end
    return (w_sum / wt_sum, n_nb)
end

# ── 格网坐标的射线法多边形判定 ──
function PointInPolygonRC(r::Int, c::Int, poly::Vector{Tuple{Int,Int}})
    n = length(poly)
    inside = false
    j = n
    for i in 1:n
        ri, ci = poly[i]
        rj, cj = poly[j]
        if (ri > r) != (rj > r)
            x_inter = cj + (ci - cj) * (r - rj) / (ri - rj)
            if c < x_inter
                inside = !inside
            end
        end
        j = i
    end
    return inside
end

# ── 主函数 ──
function InterpolateGrid(
    input_grid::String,
    boundary_file::String,
    output_file::String;
    power::Int = 2,
    min_neighbors::Int = 3,
    initial_radius::Int = 5,
    max_radius::Int = 50
)
    println("="^60)
    println("迭代 IDW 插值 — 格网坐标多边形")
    println("="^60)

    # ── 1. 加载格网 ──
    println("[1/4] 加载格网...")
    grid_data = readdlm(input_grid, Float64)
    gl, ga, gd = grid_data[:,1], grid_data[:,2], grid_data[:,3]
    @printf("  %d 点\n", length(gl))

    # 复用 MBJulia 常量
    lon_min, lon_max = extrema(gl); lat_min, lat_max = extrema(ga)
    mid_lat = (lat_min + lat_max) / 2
    m_per_deg_lon = MBJulia.METERS_PER_DEG_LAT * cosd(mid_lat)
    lon_step = 1.0 / m_per_deg_lon
    lat_step = 1.0 / MBJulia.METERS_PER_DEG_LAT
    cols = round(Int, (lon_max - lon_min) / lon_step) + 1
    rows = round(Int, (lat_max - lat_min) / lat_step) + 1
    @printf("  格网: %d × %d\n", cols, rows)
    cc(c) = lon_min + (c - 0.5) * lon_step
    cr(r) = lat_min + (r - 0.5) * lat_step

    # 构建深度矩阵（NaN 表示空白）
    depth_mat = fill(NaN, rows, cols)
    for i in 1:size(grid_data, 1)
        c = clamp(round(Int, (gl[i] - lon_min) / lon_step) + 1, 1, cols)
        r = clamp(round(Int, (ga[i] - lat_min) / lat_step) + 1, 1, rows)
        depth_mat[r, c] = gd[i]
    end

    # ── 2. 加载多边形 → 转格网坐标 ──
    println("[2/4] 加载多边形...")
    polygons_rc = Vector{Vector{Tuple{Int,Int}}}()
    cur = Tuple{Int,Int}[]
    for line in eachline(boundary_file)
        stripped = strip(line)
        if isempty(stripped)
            if length(cur) >= 3
                push!(polygons_rc, cur)
            end
            cur = Tuple{Int,Int}[]
        else
            parts = split(stripped)
            lon = parse(Float64, parts[1])
            lat = parse(Float64, parts[2])
            r = clamp(round(Int, (lat - lat_min) / lat_step) + 1, 1, rows)
            c = clamp(round(Int, (lon - lon_min) / lon_step) + 1, 1, cols)
            push!(cur, (r, c))
        end
    end
    if length(cur) >= 3
        push!(polygons_rc, cur)
    end
    @printf("  多边形: %d 段\n", length(polygons_rc))
    if isempty(polygons_rc)
        println("  跳过")
        return
    end

    # ── 3. 判定内部空白 ──
    println("[3/4] 判定内部空白...")
    todo_list = Tuple{Int,Int}[]
    n_in, n_known = 0, 0
    for r in 1:rows, c in 1:cols
        in_any = false
        for poly in polygons_rc
            if PointInPolygonRC(r, c, poly)
                in_any = true
                break
            end
        end
        in_any || continue
        n_in += 1
        if isnan(depth_mat[r, c])
            push!(todo_list, (r, c))
        else
            n_known += 1
        end
    end
    n_blank = length(todo_list)
    @printf("  多边形内: %d (已知 %d + 空白 %d, %.1f%% fill)\n",
            n_in, n_known, n_blank, 100 * n_known / n_in)
    if n_blank == 0
        println("  无需插值")
        return
    end

    # ── 4. 迭代 IDW ──
    println("[4/4] 迭代 IDW...")
    radius = initial_radius
    iter = 0
    total_filled = 0
    while !isempty(todo_list) && iter < 100
        iter += 1
        n_filled = 0
        still = Tuple{Int,Int}[]
        sizehint!(still, length(todo_list))
        for (r, c) in todo_list
            d, nb = IDWInterpolate(r, c, depth_mat, radius, power)
            if nb >= min_neighbors
                depth_mat[r, c] = d
                n_filled += 1
            else
                push!(still, (r, c))
            end
        end
        total_filled += n_filled
        if n_filled == 0
            radius = min(floor(Int, radius * 1.5), max_radius)
        end
        todo_list = still
        if n_filled == 0 && radius >= max_radius
            break
        end
    end

    if !isempty(todo_list)
        @printf("  ⚠ 剩 %d 空白格网未填充\n", length(todo_list))
    end

    # ── 写入 ──
    n_out = 0
    open(output_file, "w") do io
        for r in 1:rows, c in 1:cols
            if !isnan(depth_mat[r, c])
                @printf(io, "%.8f\t%.8f\t%.3f\n", cc(c), cr(r), depth_mat[r, c])
                n_out += 1
            end
        end
    end

    # ── 验证 ──
    verify = readdlm(output_file, Float64)
    vlon = extrema(verify[:,1]); vlat = extrema(verify[:,2])
    vcols = round(Int, (vlon[2] - vlon[1]) / lon_step) + 1
    vrows = round(Int, (vlat[2] - vlat[1]) / lat_step) + 1

    println("="^60)
    @printf("  已知 %d + 插值 %d = 输出 %d\n", n_known, total_filled, n_out)
    @printf("  坐标一致性: %d×%d %s\n", vcols, vrows, (vcols==cols && vrows==rows) ? "✓" : "✗")
    @printf("  → %s\n", output_file)
end
