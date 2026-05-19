# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Stage 5 — 第二阶段：迭代空间一致性校验                                       ║
# ║                                                                              ║
# ║  前提假设：海底地形连续 → 相邻 1m 格网水深应接近。                              ║
# ║                                                                              ║
# ║  对每个可疑格网：                                                              ║
# ║    1. 查 3×3 邻域内已解决格网的水深 → 取均值作参考值                           ║
# ║    2. 选最接近参考值的深度簇 → 标记为已解决                                    ║
# ║    3. 该格网加入已解决池 → 成为下一轮邻域参考                                  ║
# ║                                                                              ║
# ║  多轮迭代：每轮解决一批 → 邻域扩大 → 之前不够邻格的现在够了。                     ║
# ║  兜底：一轮都没解决 → 无法再推进 → 剩余全部取最大簇。                           ║
# ║                                                                              ║
# ║  优化：使用 QuestionableCell 具名结构体，消除 Dict{Any} 类型不稳定性。          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ── 收集某个格网邻域中已解决格网的参考水深 ──
function CollectNeighborDepths(key, resolved_depth::Dict, kernel)
    neighbor_depths = Float64[]
    for (dx, dy) in kernel
        nkey = (key[1] + dx, key[2] + dy)
        if haskey(resolved_depth, nkey)
            push!(neighbor_depths, resolved_depth[nkey])
        end
    end
    return neighbor_depths
end


# ── 在多个簇中选最接近参考值的一个 ──
function SelectBestCluster(clusters::Vector{UnitRange{Int}}, d_sorted::Vector{Float64},
                           ref_depth::Float64)
    best_idx = 1
    best_dist = Inf
    for (ci, cr) in enumerate(clusters)
        dist = abs(mean(d_sorted[cr]) - ref_depth)
        if dist < best_dist
            best_dist = dist
            best_idx = ci
        end
    end
    return best_idx
end


# ── 第二阶段主函数：迭代空间一致性 ──
function Stage2SpatialFilter(questionable_keys, cell_depth,
                              cell_simple, cell_questionable, n_total::Int)
    println("[5/6] 第二阶段 — 迭代空间一致性（需要 ≥ $(NEIGHBOR_K) 个已解决邻格）...")

    keep_mask = falses(n_total)
    rejected  = 0

    # 已解决池：初始 = 所有简单格网（类型安全的 Dict）
    resolved_depth = Dict{Tuple{Int,Int}, Float64}()
    for (key, idxs) in cell_simple
        resolved_depth[key] = cell_depth[key]
        keep_mask[idxs] .= true
    end

    unresolved = copy(questionable_keys)
    iteration  = 0

    while !isempty(unresolved)
        iteration += 1
        resolved_this_round = 0
        fallback_this_round = 0
        still_unresolved = Tuple{Int,Int}[]
        sizehint!(still_unresolved, length(unresolved))

        for key in unresolved
            qc = cell_questionable[key]  # QuestionableCell — 类型稳定

            neighbor_depths = CollectNeighborDepths(key, resolved_depth, KERNEL)

            if length(neighbor_depths) >= NEIGHBOR_K
                ref_depth = mean(neighbor_depths)
                best_idx = SelectBestCluster(qc.clusters, qc.d_sorted, ref_depth)
                keep_range = qc.clusters[best_idx]

                keep_mask[qc.idxs[qc.perm[keep_range]]] .= true
                rejected += length(qc.idxs) - length(keep_range)
                resolved_depth[key] = mean(qc.d_sorted[keep_range])
                resolved_this_round += 1
            else
                push!(still_unresolved, key)
            end
        end

        if resolved_this_round == 0
            # 兜底：本轮零进展 → 无法推进，剩余全部取最大簇
            for key in still_unresolved
                qc = cell_questionable[key]
                sizes = [length(r) for r in qc.clusters]
                keep_range = qc.clusters[argmax(sizes)]
                keep_mask[qc.idxs[qc.perm[keep_range]]] .= true
                rejected += length(qc.idxs) - length(keep_range)
                fallback_this_round += 1
            end
            @printf("  第 %d 轮: 解决 %d，兜底 %d（无法继续推进 → 退出）\n",
                    iteration, resolved_this_round, fallback_this_round)
            break
        end

        @printf("  第 %d 轮: 解决 %d，仍待定 %d\n",
                iteration, resolved_this_round, length(still_unresolved))
        unresolved = still_unresolved
    end

    n_kept = count(keep_mask)
    @printf("  收敛于 %d 轮。保留: %d / %d (%.1f%%)，剔除: %d\n",
            iteration, n_kept, n_total, n_kept / n_total * 100,
            n_total - n_kept)

    return keep_mask, rejected
end
