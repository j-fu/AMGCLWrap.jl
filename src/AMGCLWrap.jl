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

struct xxxSolver{Tv, Ti}
    handle::Ptr{Cvoid}
    blocksize::Cint
end


function xxxtest(x::Tv, n::Ti) where{Tv,Ti}
    s=ccall( (:xxxCreate, libamgcl_c), xxxSolver{Tv,Ti}, (Cint,), n)
    s.blocksize
end

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon

end

