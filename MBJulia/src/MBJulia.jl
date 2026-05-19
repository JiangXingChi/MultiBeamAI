# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  MBJulia — 多波束水深数据两阶段迭代滤波                                       ║
# ║                                                                              ║
# ║  模块入口。include 各子文件，仅导出 Run 函数。                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

module MBJulia

using DelimitedFiles
using Printf
using Statistics

# ── 常量 (Const) ──
include("01_constants.jl")

# ── 数据结构 (Types) ──
include("02_types.jl")

# ── 各处理阶段 ──
include("03_io.jl")         # Stage 1 (Load) + Stage 6 (Write)
include("04_geometry.jl")   # Stage 2
include("05_binning.jl")    # Stage 3
include("06_stage1.jl")     # Stage 4
include("07_stage2.jl")     # Stage 5
include("08_pipeline.jl")   # Orchestrator (一次处理)

# 仅暴露编排入口和共享工具
export Run, ComputeGridGeometry

end
