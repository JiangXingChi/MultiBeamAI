# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  边界检测 — 迭代剥离法                                                           ║
# ║                                                                              ║
# ║  算法：                                                                        ║
# ║    1. 连通区域标记 (8-邻域 BFS)                                               ║
# ║    2. 取最大连通区域 → 边界格网检测 → 边界追踪 → DP简化 → 输出多边形              ║
# ║    3. 射线法判定多边形内部 → 删除内部所有点                                      ║
# ║    4. 回到步骤1，直到无剩余连通区域                                              ║
# ║                                                                              ║
# ║  优化：Moore-neighbor 追踪替代 O(n²) 半径搜索；缓存边界计数；                     ║
# ║         复用 MBJulia 常量而非重复声明。                                         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using DelimitedFiles
using Printf

const DIRS8 = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

# 8个方向的 Moore-neighbor 搜索顺序（顺时针，从上方开始）
const MOORE_DIRS = [(-1,0), (-1,1), (0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1)]

function PointInPolygon(px::Float64, py::Float64,
                         poly_x::Vector{Float64}, poly_y::Vector{Float64})::Bool
    n = length(poly_x); inside = false; j = n
    for i in 1:n
        if ((poly_y[i] > py) != (poly_y[j] > py)) &&
           (px < (poly_x[j]-poly_x[i])*(py-poly_y[i])/(poly_y[j]-poly_y[i]) + poly_x[i])
            inside = !inside
        end; j = i
    end; return inside
end

# ── Moore-neighbor 边界追踪（快速，O(n)）──
# 返回追踪到的边界点 (poly_r, poly_c)，失败返回空
function MooreTrace(boundary::BitMatrix, start_r::Int, start_c::Int,
                     rows::Int, cols::Int)
    visited = falses(rows, cols)
    visited[start_r, start_c] = true
    poly_r = Int[start_r]; poly_c = Int[start_c]

    curr_r, curr_c = start_r, start_c
    # 初始进入方向：假设从左边进入
    entry_dir = 1  # index into MOORE_DIRS

    max_steps = min(count(boundary) * 3, 100_000)  # 安全上限
    steps = 0

    while true
        steps += 1
        steps > max_steps && break

        found = false
        # 从 entry_dir 的逆时针侧开始搜索（Moore-neighbor 标准做法）
        search_start = mod(entry_dir + 5, 8) + 1  # 顺时针转约180°开始

        for offset in 0:7
            di = mod(search_start + offset - 1, 8) + 1
            dr, dc = MOORE_DIRS[di]
            nr, nc = curr_r + dr, curr_c + dc

            if 1 <= nr <= rows && 1 <= nc <= cols && boundary[nr, nc] && !visited[nr, nc]
                visited[nr, nc] = true
                push!(poly_r, nr); push!(poly_c, nc)
                curr_r, curr_c = nr, nc
                entry_dir = di
                found = true
                break
            end
        end

        if !found
            break  # 追踪结束（回到起点或死胡同）
        end

        # 检查是否回到起点附近（闭合）
        if length(poly_r) > 8
            dr_s = abs(curr_r - start_r)
            dc_s = abs(curr_c - start_c)
            if dr_s <= 1 && dc_s <= 1
                push!(poly_r, start_r); push!(poly_c, start_c)
                break
            end
        end
    end

    return (poly_r, poly_c)
end

# ── 备用：最近邻追踪（当 Moore 追踪失败时）──
function NearestTrace(boundary::BitMatrix, start_r::Int, start_c::Int,
                       rows::Int, cols::Int, mid_lat::Float64, n_boundary::Int)
    visited = falses(rows, cols)
    visited[start_r, start_c] = true
    poly_r = Int[start_r]; poly_c = Int[start_c]
    curr_r, curr_c = start_r, start_c

    cos_fac = cosd(mid_lat)^2
    max_steps = min(n_boundary * 2, 50_000)

    for _ in 1:max_steps
        best_d2 = Inf; best_r = best_c = 0
        # 搜索所有边界格网中最近的未访问邻格
        for r in max(1, curr_r-50):min(rows, curr_r+50)
            for c in max(1, curr_c-50):min(cols, curr_c+50)
                boundary[r, c] || continue
                visited[r, c] && continue
                d2 = (r - curr_r)^2 + (c - curr_c)^2 * cos_fac
                if d2 < best_d2
                    best_d2 = d2; best_r, best_c = r, c
                end
            end
        end
        if best_d2 == Inf || best_d2 > 2500  # > 50格网距离
            break
        end
        curr_r, curr_c = best_r, best_c
        visited[curr_r, curr_c] = true
        push!(poly_r, curr_r); push!(poly_c, curr_c)
    end

    push!(poly_r, start_r); push!(poly_c, start_c)
    return (poly_r, poly_c)
end

# ── Douglas-Peucker 简化 ──
function DPSimplify(xs, ys, i_start, i_end, epsilon)
    dmax = 0.0; imax = i_start
    if i_end > i_start + 1
        x0, y0 = xs[i_start], ys[i_start]; x1, y1 = xs[i_end], ys[i_end]
        seg2 = (x1 - x0)^2 + (y1 - y0)^2
        for k in i_start+1:i_end-1
            d = if seg2 < 1e-20
                sqrt((xs[k] - x0)^2 + (ys[k] - y0)^2)
            else
                t = clamp(((xs[k]-x0)*(x1-x0)+(ys[k]-y0)*(y1-y0))/seg2, 0.0, 1.0)
                sqrt((xs[k]-(x0+t*(x1-x0)))^2 + (ys[k]-(y0+t*(y1-y0)))^2)
            end
            if d > dmax; dmax = d; imax = k; end
        end
    end
    if dmax > epsilon
        left  = DPSimplify(xs, ys, i_start, imax, epsilon)
        right = DPSimplify(xs, ys, imax, i_end, epsilon)
        return [left[1:end-1]; right]
    else
        return [i_start, i_end]
    end
end

# ── 主函数 ──
function DetectBoundary(input_file::String, output_file::String;
                         min_region_cells::Int = 10)
    println("="^60)
    println("边界检测 — 迭代剥离法")
    println("="^60)

    # ── 加载 ──
    raw = readdlm(input_file, Float64; comments=false)
    all_lons = raw[:,1]; all_lats = raw[:,2]; all_depths = raw[:,3]
    n_total = length(all_lons)
    @printf("  %d 点, 深度 %.2f ~ %.2f m\n", n_total,
            minimum(all_depths), maximum(all_depths))

    # 格网几何 — 复用 MBJulia 常量
    lon_min, lon_max = extrema(all_lons); lat_min, lat_max = extrema(all_lats)
    mid_lat = (lat_min + lat_max) / 2
    m_per_deg_lon = MBJulia.METERS_PER_DEG_LAT * cosd(mid_lat)
    lon_step = 1.0 / m_per_deg_lon
    lat_step = 1.0 / MBJulia.METERS_PER_DEG_LAT
    cols = round(Int, (lon_max - lon_min) / lon_step) + 1
    rows = round(Int, (lat_max - lat_min) / lat_step) + 1
    @printf("  格网: %d × %d\n", cols, rows)
    cc(c) = lon_min + (c - 0.5) * lon_step
    cr(r) = lat_min + (r - 0.5) * lat_step

    # 活跃点掩码
    active = trues(n_total)

    all_polygons = Vector{Float64}[]
    total_vertices = 0
    iteration = 0
    max_iterations = 20

    EPSILON_DEG = 1.0 / MBJulia.METERS_PER_DEG_LAT

    while true
        iteration += 1
        if iteration > max_iterations
            break
        end

        # ── 构建当前活跃点的掩膜 ──
        mask = falses(rows, cols)
        for i in 1:n_total
            if active[i]
                c = clamp(round(Int, (all_lons[i]-lon_min)/lon_step) + 1, 1, cols)
                r = clamp(round(Int, (all_lats[i]-lat_min)/lat_step) + 1, 1, rows)
                mask[r, c] = true
            end
        end

        # ── 连通区域标记 ──
        labels = zeros(Int, rows, cols)
        region_sizes = Int[]
        next_label = 0
        for sr in 1:rows, sc in 1:cols
            if mask[sr, sc] && labels[sr, sc] == 0
                next_label += 1
                queue = [(sr, sc)]
                labels[sr, sc] = next_label
                cnt = 1
                while !isempty(queue)
                    r, c = popfirst!(queue)
                    for (dr, dc) in DIRS8
                        nr, nc = r + dr, c + dc
                        if 1 <= nr <= rows && 1 <= nc <= cols && mask[nr, nc] && labels[nr, nc] == 0
                            labels[nr, nc] = next_label
                            push!(queue, (nr, nc))
                            cnt += 1
                        end
                    end
                end
                push!(region_sizes, cnt)
            end
        end

        valid = findall(s -> s >= min_region_cells, region_sizes)
        isempty(valid) && break

        # 取最大区域
        max_idx = valid[argmax(region_sizes[valid])]
        max_cells = region_sizes[max_idx]
        @printf("\n  迭代 %d: 最大区域 #%d (%d cells)\n", iteration, max_idx, max_cells)

        # ── 对该区域做边界检测 ──
        region_mask = labels .== max_idx
        boundary = falses(rows, cols)
        n_boundary = 0
        for r in 1:rows, c in 1:cols
            region_mask[r, c] || continue
            is_b = false
            for (dr, dc) in DIRS8
                nr, nc = r + dr, c + dc
                if nr < 1 || nr > rows || nc < 1 || nc > cols || !region_mask[nr, nc]
                    is_b = true; break
                end
            end
            if is_b
                boundary[r, c] = true
                n_boundary += 1
            end
        end

        @printf("    边界格网: %d\n", n_boundary)

        # 寻找起始点
        start_r = start_c = 0
        for c2 in 1:cols
            for r2 in 1:rows
                if boundary[r2, c2]
                    start_r, start_c = r2, c2
                    break
                end
            end
            start_r != 0 && break
        end
        if start_r == 0
            start_r, start_c = findfirst(boundary)
        end

        # ── 追踪（优先 Moore-neighbor，失败回退最近邻）──
        poly_r, poly_c = MooreTrace(boundary, start_r, start_c, rows, cols)
        if length(poly_r) < 4
            poly_r, poly_c = NearestTrace(boundary, start_r, start_c, rows, cols, mid_lat, n_boundary)
        end

        # 转换为坐标
        xs = [cc(poly_c[i]) for i in 1:length(poly_r)]
        ys = [cr(poly_r[i]) for i in 1:length(poly_r)]

        # DP 简化
        idx = DPSimplify(xs, ys, 1, length(xs), EPSILON_DEG)
        poly_x = [xs[idx[i]] for i in 1:length(idx)]
        poly_y = [ys[idx[i]] for i in 1:length(idx)]

        @printf("    追踪 %d → DP %d 顶点\n", length(poly_r), length(idx))

        # 构建多边形向量
        poly = vcat(poly_x, poly_y)
        push!(all_polygons, poly)
        total_vertices += length(idx)

        # ── 剥离：删除多边形内部的所有活跃点 ──
        n_removed = 0
        for i in 1:n_total
            if active[i] && PointInPolygon(all_lons[i], all_lats[i], poly_x, poly_y)
                active[i] = false
                n_removed += 1
            end
        end
        @printf("    剥离 %d 点, 剩余 %d\n", n_removed, count(active))

        if n_removed == 0
            println("    多边形未剥离任何点，终止迭代")
            pop!(all_polygons)
            break
        end
    end

    # ── 写入 ──
    println("\n写入...")
    open(output_file, "w") do io
        for (i, p) in enumerate(all_polygons)
            n = length(p) ÷ 2
            for j in 1:n
                @printf(io, "%.8f\t%.8f\n", p[j], p[n+j])
            end
            i < length(all_polygons) && println(io)
        end
    end

    println("="^60)
    println("完成")
    @printf("  多边形: %d 段, 总顶点: %d\n", length(all_polygons), total_vertices)
    @printf("  → %s\n", output_file)
end
