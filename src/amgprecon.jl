mutable struct AMGPrecon{Tv,Ti}
    handle::Ptr{Cvoid}
end

function AMGPrecon(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLAMGPreconCreate",libamgcl_c),
               AMGPrecon{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cstring),
               n, ia,ja, a,param);
    finalizer(s->ccall(("amgclcDLAMGPreconDestroy",libamgcl_c),
                       Cvoid,
                       (AMGPrecon{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::AMGPrecon{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLAMGPreconApply",libamgcl_c),
          AMGCLInfo,
          (AMGPrecon{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end

LinearAlgebra.ldiv!(s::AMGPrecon, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::AMGPrecon, v)
    apply!(s,u,v)
    u
end

    
function AMGPrecon(csr::SparseMatrixCSR{Bi,Tv,Ti},param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=AMGPrecon(csr.m, csr.rowptr,csr.colval,csr.nzval,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function AMGPrecon(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    AMGPrecon(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),param)
end

