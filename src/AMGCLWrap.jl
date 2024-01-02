module AMGCLWrap
using AMGCL_C_jll
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
import JSON3

include("common.jl")
include("amgsolver.jl")
include("rlxsolver.jl")
include("amgprecon.jl")
include("rlxprecon.jl")

include("blockamgsolver.jl")
include("blockrlxsolver.jl")
include("blockamgprecon.jl")
include("blockrlxprecon.jl")

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon
export BlockAMGSolver, BlockRLXSolver, BlockAMGPrecon, BlockRLXPrecon

end
