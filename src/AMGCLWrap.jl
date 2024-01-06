module AMGCLWrap
using AMGCL_C_jll
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
import JSON3


include("amgclc_wrapper.jl")

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon

end

