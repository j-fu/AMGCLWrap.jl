mutable struct BlockRLXPrecon{Tv,Ti}
    handle::Ptr{Cvoid}
    blocksize::Cint
end

function BlockRLXPrecon(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv},blocksize, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLBlockRLXPreconCreate",libamgcl_c),
               BlockRLXPrecon{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cint,Cstring),
               n, ia,ja, a, bkocksize, param);
    finalizer(s->ccall(("amgclcDLBlockRLXPreconDestroy",libamgcl_c),
                       Cvoid,
                       (BlockRLXPrecon{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::BlockRLXPrecon{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLBlockRLXPreconApply",libamgcl_c),
          AMGCLInfo,
          (BlockRLXPrecon{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end

LinearAlgebra.ldiv!(s::BlockRLXPrecon, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::BlockRLXPrecon, v)
    apply!(s,u,v)
    u
end

    
function BlockRLXPrecon(csr::SparseMatrixCSR{Bi,Tv,Ti},blocksize,param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=BlockRLXPrecon(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function BlockRLXPrecon(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},blocksize,param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    BlockRLXPrecon(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),blocksize,param)
end

