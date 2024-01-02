mutable struct BlockAMGSolver{Tv,Ti}
    handle::Ptr{Cvoid}
    blocksize::Cint
end

function BlockAMGSolver(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, blocksize, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLBlockAMGSolverCreate",libamgcl_c),
               BlockAMGSolver{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cint, Cstring),
               n, ia,ja, a, blocksize, param);
    finalizer(s->ccall(("amgclcDLBlockAMGSolverDestroy",libamgcl_c),
                       Cvoid,
                       (BlockAMGSolver{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::BlockAMGSolver{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLBlockAMGSolverApply",libamgcl_c),
          AMGCLInfo,
          (BlockAMGSolver{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end

LinearAlgebra.ldiv!(s::BlockAMGSolver, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::BlockAMGSolver, v)
    apply!(s,u,v)
    u
end
LinearAlgebra.:\(s::BlockAMGSolver, v) = ldiv!(copy(v),s,v)
    
function BlockAMGSolver(csr::SparseMatrixCSR{Bi,Tv,Ti},blocksize,param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=BlockAMGSolver(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function BlockAMGSolver(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},blocksize,param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    BlockAMGSolver(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),blocksize,param)
end

