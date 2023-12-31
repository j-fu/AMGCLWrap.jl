mutable struct AMGSolver{Tv,Ti}
    handle::Ptr{Cvoid}
end

function AMGSolver(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLAMGSolverCreate",libamgcl_c),
               AMGSolver{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cstring),
               n, ia,ja, a,param);
    finalizer(s->ccall(("amgclcDLAMGSolverDestroy",libamgcl_c),
                       Cvoid,
                       (AMGSolver{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::AMGSolver{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLAMGSolverApply",libamgcl_c),
          AMGCLInfo,
          (AMGSolver{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end

LinearAlgebra.ldiv!(s::AMGSolver, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::AMGSolver, v)
    apply!(s,u,v)
    u
end
LinearAlgebra.:\(s::AMGSolver, v) = ldiv!(copy(v),s,v)
    
function AMGSolver(csr::SparseMatrixCSR{Bi,Tv,Ti},param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=AMGSolver(csr.m, csr.rowptr,csr.colval,csr.nzval,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function AMGSolver(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    AMGSolver(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),param)
end

