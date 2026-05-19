# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  数据结构 — 格网几何                                                          ║
# ║                                                                              ║
# ║  为什么用 struct：格网的度数步长、行列数、边界等 12 个值彼此关联，              ║
# ║  散落为裸变量容易在传递时遗漏或混淆。封装后一个变量承载全部几何信息。            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

struct GridGeometry
    lon_min::Float64
    lon_max::Float64
    lat_min::Float64
    lat_max::Float64
    mid_lat::Float64                   # 测区中心纬度 → 用于 cos 校正
    meters_per_deg_lon::Float64        # 该纬度处 1° 经度的地面距离
    lon_step::Float64                  # 1m 经度步长（度）
    lat_step::Float64                  # 1m 纬度步长（度）
    width_m::Float64                   # 东西跨度（米）
    height_m::Float64                  # 南北跨度（米）
    cols::Int                          # 格网列数
    rows::Int                          # 格网行数
end

# 可疑格网数据 — 第二阶段投票所需信息
# 用具名结构体替代 NamedTuple/Any，消除类型不稳定性
struct QuestionableCell
    idxs::Vector{Int}                  # 该格网内所有点的全局索引
    perm::Vector{Int}                  # 排序排列 (sortperm)
    clusters::Vector{UnitRange{Int}}   # 深度簇区间列表
    d_sorted::Vector{Float64}          # 排序后的深度值
end
