mutable struct RLXPrecon{Tv,Ti}
    handle::Ptr{Cvoid}
end

function RLXPrecon(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, param::String) where {Tv<:Float64,Ti<:Int64}
    this=ccall(("amgclcDLRLXPreconCreate",libamgcl_c),
               RLXPrecon{Float64,Int64},
               (Cint, Ptr{Clong}, Ptr{Clong},Ptr{Cdouble},Cstring),
               n, ia,ja, a,param);
    finalizer(s->ccall(("amgclcDLRLXPreconDestroy",libamgcl_c),
                       Cvoid,
                       (RLXPrecon{Float64, Int64},),
                       s),this)
    this
end

function apply!(s::RLXPrecon{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:Float64,Ti<:Int64}
    i=ccall(("amgclcDLRLXPreconApply",libamgcl_c),
          AMGCLInfo,
          (RLXPrecon{Float64, Int64},Ptr{Cdouble}, Ptr{Cdouble}),
          s,sol,rhs)
    (iters=i.iters, residual=i.residual)
end


LinearAlgebra.ldiv!(s::RLXPrecon, v) = v.=ldiv!(copy(v),s,v)
function LinearAlgebra.ldiv!(u, s::RLXPrecon, v)
    apply!(s,u,v)
    u
end




function RLXPrecon(csr::SparseMatrixCSR{Bi,Tv,Ti},param=nothing) where {Bi,Tv,Ti}
    if csr.m!=csr.n
        error("Matrix must be square")
    end
    myoffset=1-getoffset(Bi)
    csr.rowptr.-=myoffset
    csr.colval.-=myoffset
    s=RLXPrecon(csr.m, csr.rowptr,csr.colval,csr.nzval,myjson(param))
    csr.rowptr.+=myoffset
    csr.colval.+=myoffset
    return s
end

function RLXPrecon(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},param="") where {Tv,Ti}
    if csc.m!=csc.n
        error("Matrix must be square")
    end
    RLXPrecon(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),param)
end
