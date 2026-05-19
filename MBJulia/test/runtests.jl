# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  测试占位 — 功能测试后续添加                                                    ║
# ║                                                                              ║
# ║  用法:                                                                        ║
# ║    julia --project=. -e 'using Pkg; Pkg.test()'                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using MBJulia
using Test

@testset "Constants" begin
    @test MBJulia.DEPTH_GAP == 0.25
    @test MBJulia.NEIGHBOR_K == 3
    @test length(MBJulia.KERNEL) == 8
    @test MBJulia.METERS_PER_DEG_LAT ≈ 111_320.0
end

@testset "SplitByDepthGaps" begin
    # 无间隙：等间距数据
    d = collect(1.0:0.1:10.0)
    @test isempty(MBJulia.SplitByDepthGaps(d, 1.0))

    # 有间隙：中间跳变
    d2 = [1.0, 1.1, 5.0, 5.1, 10.0, 10.1]
    gaps = MBJulia.SplitByDepthGaps(d2, 1.0)
    @test gaps == [3, 5]
end

@testset "BuildClustersFromGaps" begin
    clusters = MBJulia.BuildClustersFromGaps([3, 5], 6)
    @test length(clusters) == 3
    @test clusters[1] == 1:2
    @test clusters[2] == 3:4
    @test clusters[3] == 5:6
end

@testset "SelectBestCluster" begin
    d_sorted = [1.0, 1.1, 3.0, 3.1, 9.0, 9.1]
    clusters = [1:2, 3:4, 5:6]   # 均值 ~1.05, ~3.05, ~9.05
    @test MBJulia.SelectBestCluster(clusters, d_sorted, 3.0) == 2   # 3.05 最近 3.0
    @test MBJulia.SelectBestCluster(clusters, d_sorted, 9.0) == 3   # 9.05 最近 9.0
    @test MBJulia.SelectBestCluster(clusters, d_sorted, 1.0) == 1   # 1.05 最近 1.0
end

println("All tests passed.")
