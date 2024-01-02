mutable struct BlockAMGPrecon{Tv,Ti}
    handle::Ptr{Cvoid}
    blocksize::Cint
end

function BlockAMGPrecon(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv},blocksize, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLBlockAMGPreconCreate",libamgcl_c),
               BlockAMGPrecon{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cint,Cstring),
               n, ia,ja, a, bkocksize, param);
    finalizer(s->ccall(("amgclcDLBlockAMGPreconDestroy",libamgcl_c),
                       Cvoid,
                       (BlockAMGPrecon{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::BlockAMGPrecon{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLBlockAMGPreconApply",libamgcl_c),
          AMGCLInfo,
          (BlockAMGPrecon{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end

LinearAlgebra.ldiv!(s::BlockAMGPrecon, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::BlockAMGPrecon, v)
    apply!(s,u,v)
    u
end

    
function BlockAMGPrecon(csr::SparseMatrixCSR{Bi,Tv,Ti},blocksize,param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=BlockAMGPrecon(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function BlockAMGPrecon(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},blocksize,param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    BlockAMGPrecon(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),blocksize,param)
end

