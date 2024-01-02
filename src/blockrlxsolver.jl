mutable struct BlockRLXSolver{Tv,Ti}
    handle::Ptr{Cvoid}
    blocksize::Cint
end

function BlockRLXSolver(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, blocksize, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLBlockRLXSolverCreate",libamgcl_c),
               BlockRLXSolver{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cint,Cstring),
               n, ia,ja, a, blocksize, param);
    finalizer(s->ccall(("amgclcDLBlockRLXSolverDestroy",libamgcl_c),
                       Cvoid,
                       (BlockRLXSolver{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::BlockRLXSolver{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLBlockRLXSolverApply",libamgcl_c),
          AMGCLInfo,
          (BlockRLXSolver{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end


LinearAlgebra.ldiv!(s::BlockRLXSolver, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::BlockRLXSolver, v)
    apply!(s,u,v)
    u
end
LinearAlgebra.:\(s::BlockRLXSolver, v) = ldiv!(copy(v),s,v)



function BlockRLXSolver(csr::SparseMatrixCSR{Bi,Tv,Ti},blocksize, param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=BlockRLXSolver(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function BlockRLXSolver(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},blocksize, param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    BlockRLXSolver(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),blocksize, param)
end
