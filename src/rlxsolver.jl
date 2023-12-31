mutable struct RLXSolver{Tv,Ti}
    handle::Ptr{Cvoid}
end

function RLXSolver(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLRLXSolverCreate",libamgcl_c),
               RLXSolver{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cstring),
               n, ia,ja, a,param);
    finalizer(s->ccall(("amgclcDLRLXSolverDestroy",libamgcl_c),
                       Cvoid,
                       (RLXSolver{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::RLXSolver{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLRLXSolverApply",libamgcl_c),
          AMGCLInfo,
          (RLXSolver{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end


LinearAlgebra.ldiv!(s::RLXSolver, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::RLXSolver, v)
    apply!(s,u,v)
    u
end
LinearAlgebra.:\(s::RLXSolver, v) = ldiv!(copy(v),s,v)



function RLXSolver(csr::SparseMatrixCSR{Bi,Tv,Ti},param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=RLXSolver(csr.m, csr.rowptr,csr.colval,csr.nzval,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function RLXSolver(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    RLXSolver(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),param)
end
