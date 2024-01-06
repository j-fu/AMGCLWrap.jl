module AMGCLWrap
using AMGCL_C_jll
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
import JSON3


include("amgclc_wrapper.jl")
function simpletest()
    ccall( (:simpletest, libamgcl_c), Cint, ())
end

function fulltest(n)
    ccall( (:fulltest, libamgcl_c), Cint, (Cint,), n)
end

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon

end

