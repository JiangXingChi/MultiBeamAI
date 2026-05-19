# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Stage 4 — 第一阶段：格网内深度间隙聚类                                       ║
# ║                                                                              ║
# ║  核心思路：对每个格网内的深度值排序后扫描相邻间隙。                              ║
# ║  如果深度跨度 < DEPTH_GAP → 一种水深 → 优质（直接取均值）                       ║
# ║  如果有间隙 → 在间隙处切分成多个深度簇：                                        ║
# ║    · 恰好 2 簇 + 大簇 > 50% → 海底 + 少量噪声 → 优质                          ║
# ║    · ≥3 簇 或 两簇大小接近 → 无法自动判断 → 可疑 → 移交第二阶段                ║
# ║                                                                              ║
# ║  本质：一维 DBSCAN 的简化版（ε = DEPTH_GAP），不做 MinPts 检查。                ║
# ║                                                                              ║
# ║  优化：用具名 Struct 替代异构 NamedTuple，消除 Dict{Any} 类型不稳定性。          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ── 纯函数：在排序深度序列中找出所有 > gap 的断点 ──
function SplitByDepthGaps(d_sorted::Vector{Float64}, gap::Float64)
    gaps = Int[]
    for j in 2:length(d_sorted)
        if d_sorted[j] - d_sorted[j-1] > gap
            push!(gaps, j)
        end
    end
    return gaps
end


# ── 纯函数：根据间隙位置将深度序列切分为多个簇 ──
function BuildClustersFromGaps(gaps::Vector{Int}, n_total::Int)
    clusters = Vector{UnitRange{Int}}()
    sizehint!(clusters, length(gaps) + 1)
    start_idx = 1
    for g in gaps
        push!(clusters, start_idx:(g-1))
        start_idx = g
    end
    push!(clusters, start_idx:n_total)
    return clusters
end


# ── 对一个格网做分簇 ──
# 返回 (status::Symbol, depth::Float64, questionable::Union{QuestionableCell,Nothing})
# status = :simple  → 可直接用 depth 值
# status = :questionable → 需要第二阶段空间投票，信息在 questionable 中
function ClassifyCell(d::Vector{Float64}, idxs::Vector{Int}, gap::Float64)
    d_range = maximum(d) - minimum(d)

    if d_range < gap
        return (:simple, mean(d), nothing)
    end

    perm = sortperm(d)
    d_sorted = d[perm]
    gaps_found = SplitByDepthGaps(d_sorted, gap)

    if isempty(gaps_found)
        return (:simple, mean(d), nothing)
    end

    clusters = BuildClustersFromGaps(gaps_found, length(d_sorted))
    qcell = QuestionableCell(idxs, perm, clusters, d_sorted)
    return (:questionable, NaN, qcell)
end


# ── 第一阶段主函数：对所有格网执行间隙聚类 ──
function Stage1GapFilter(cell_to_indices, depths)
    println("[4/6] 第一阶段 — 格网内分簇（间隙阈值 = $(DEPTH_GAP)m）...")

    # 分离存储：简单格网 vs 可疑格网，消除 Dict{Any} 类型不稳定性
    cell_depth   = Dict{Tuple{Int,Int}, Float64}()
    cell_simple  = Dict{Tuple{Int,Int}, Vector{Int}}()
    cell_questionable = Dict{Tuple{Int,Int}, QuestionableCell}()
    questionable_keys = Tuple{Int,Int}[]

    @time for (key, idxs) in cell_to_indices
        status, depth, qcell = ClassifyCell(depths[idxs], idxs, DEPTH_GAP)

        if status == :simple
            cell_depth[key]  = depth
            cell_simple[key] = qcell === nothing ? idxs : error("unexpected")
        else
            cell_depth[key]        = NaN
            cell_questionable[key] = qcell
            push!(questionable_keys, key)
        end
    end

    n_simple = length(cell_simple)
    n_questionable = length(questionable_keys)
    @printf("  简单格网: %d，含间隙格网: %d\n", n_simple, n_questionable)
    return cell_depth, cell_simple, cell_questionable, questionable_keys
end
