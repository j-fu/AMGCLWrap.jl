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

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon

end
