module AMGCLWrap
using AMGCL_C_jll
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
import JSON3
using CompositeStructs: @composite
using DocStringExtensions

include("amgclc_wrapper.jl")
include("parameters.jl")
include("struct_api.jl")

export AMGSolver, RLXSolver, AMGPrecon, RLXPrecon, blocksize_instantiated, error_state
export AMGSolverAlgorithm, RLXSolverAlgorithm


end
