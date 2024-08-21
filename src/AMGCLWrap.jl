"""
    AMGCLWrap

$(read(joinpath(@__DIR__,"..","README.md"),String))
"""
module AMGCLWrap

using AMGCL_C_jll: AMGCL_C_jll, libamgcl_c
using DocStringExtensions: DocStringExtensions, TYPEDEF, TYPEDFIELDS
using LinearAlgebra: LinearAlgebra, issymmetric, ldiv!, transpose
using SparseArrays: SparseArrays, AbstractSparseMatrix, SparseMatrixCSC, sparse
using SparseMatricesCSR: SparseMatricesCSR, SparseMatrixCSR, getoffset
using CompositeStructs: @composite
import JSON3

include("amgclc_wrapper.jl")
include("parameters.jl")
include("struct_api.jl")

export AMGSolver, RLXSolver, AMGPrecon, RLXPrecon, blocksize_instantiated, error_state
export AMGPreconditioner, RLXPreconditioner
export AMGSolverAlgorithm, RLXSolverAlgorithm


end
