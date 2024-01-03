module AMGCLWrap
using AMGCL_C_jll
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
import JSON3

struct AMGCLInfo
    iters::Cint
    residual::Cdouble
end

tojson(param)=JSON3.write(param)
tojson(::Nothing)=""
tojson(s::String)=s

const DLAMGSolver=Dict(
    "XXXMethod" => "AMGSolver",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLRLXSolver=Dict(
    "XXXMethod" => "RLXSolver",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLAMGPrecon=Dict(
    "XXXMethod" => "AMGPrecon",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLRLXPrecon=Dict(
    "XXXMethod" => "RLXPrecon",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLBlockAMGSolver=Dict(
    "XXXMethod" => "BlockAMGSolver",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLBlockRLXSolver=Dict(
    "XXXMethod" => "BlockRLXSolver",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLBlockAMGPrecon=Dict(
    "XXXMethod" => "BlockAMGPrecon",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const DLBlockRLXPrecon=Dict(
    "XXXMethod" => "BlockRLXPrecon",
    "TvTi" => "DL",
    "JTv"=> "Float64",
    "JTi"=> "Int64",
    "CTv" => "Cdouble",
    "CTi" => "Clong")

const scalar_instances=[DLAMGSolver, DLRLXSolver, DLAMGPrecon,DLRLXPrecon]
const block_instances=[DLBlockAMGSolver, DLBlockRLXSolver, DLBlockAMGPrecon,DLBlockRLXPrecon]

for instance in scalar_instances
    XXXMethod=Symbol(instance["XXXMethod"])
    TvTiXXXMethod=instance["TvTi"]*instance["XXXMethod"]
    amgclcTvTiXXXMethodCreate="amgclc$(TvTiXXXMethod)Create"
    amgclcTvTiXXXMethodDestroy="amgclc$(TvTiXXXMethod)Destroy"
    amgclcTvTiXXXMethodApply="amgclc$(TvTiXXXMethod)Apply"
    JTv=Symbol(instance["JTv"])
    JTi=Symbol(instance["JTi"])
    CTv=Symbol(instance["CTv"])
    CTi=Symbol(instance["CTi"])

    @eval begin
        mutable struct $XXXMethod{Tv,Ti}
            handle::Ptr{Cvoid}
        end
        
        function $XXXMethod(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, param::String) where {Tv<:$JTv,Ti<:$JTi}
            this=ccall(($amgclcTvTiXXXMethodCreate,libamgcl_c),
                       $XXXMethod{$JTv,$JTi},
                       (Cint, Ptr{$CTi}, Ptr{$CTi},Ptr{$CTv},Cstring),
                       n, ia,ja, a,param);
            finalizer(s->ccall(($amgclcTvTiXXXMethodDestroy,libamgcl_c),
                               Cvoid,
                               ($XXXMethod{$JTv, $JTi},),
                               s),this)
            this
        end
        
        function apply!(s::$XXXMethod{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:$JTv,Ti<:$JTi}
            i=ccall(($amgclcTvTiXXXMethodApply,libamgcl_c),
                    AMGCLInfo,
                    ($XXXMethod{$JTv, $JTi},Ptr{$CTv}, Ptr{$CTv}),
                    s,sol,rhs)
            (iters=i.iters, residual=i.residual)
        end
    end
end

for instance in block_instances
    XXXMethod=Symbol(instance["XXXMethod"])
    TvTiXXXMethod=instance["TvTi"]*instance["XXXMethod"]
    amgclcTvTiXXXMethodCreate="amgclc$(TvTiXXXMethod)Create"
    amgclcTvTiXXXMethodDestroy="amgclc$(TvTiXXXMethod)Destroy"
    amgclcTvTiXXXMethodApply="amgclc$(TvTiXXXMethod)Apply"
    JTv=Symbol(instance["JTv"])
    JTi=Symbol(instance["JTi"])
    CTv=Symbol(instance["CTv"])
    CTi=Symbol(instance["CTi"])

    @eval begin
        mutable struct $XXXMethod{Tv,Ti}
            handle::Ptr{Cvoid}
            blocksize::Cint
        end
        
        function $XXXMethod(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, blocksize,param::String) where {Tv<:$JTv,Ti<:$JTi}
            this=ccall(($amgclcTvTiXXXMethodCreate,libamgcl_c),
                       $XXXMethod{$JTv,$JTi},
                       (Cint, Ptr{$CTi}, Ptr{$CTi},Ptr{$CTv},Cint,Cstring),
                       n, ia,ja, a, blocksize, param);
            finalizer(s->ccall(($amgclcTvTiXXXMethodDestroy,libamgcl_c),
                               Cvoid,
                               ($XXXMethod{$JTv, $JTi},),
                               s),this)
            this
        end
        
        function apply!(s::$XXXMethod{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:$JTv,Ti<:$JTi}
            i=ccall(($amgclcTvTiXXXMethodApply,libamgcl_c),
                    AMGCLInfo,
                    ($XXXMethod{$JTv, $JTi},Ptr{$CTv}, Ptr{$CTv}),
                    s,sol,rhs)
            (iters=i.iters, residual=i.residual)
        end
    end
end


for XXXMethod in unique([Symbol(instance["XXXMethod"]) for instance in scalar_instances])
    @eval begin
        LinearAlgebra.ldiv!(s::$XXXMethod, v) = v.=ldiv!(copy(v),s,v)
        function LinearAlgebra.ldiv!(u, s::$XXXMethod, v)
            apply!(s,u,v)
            u
        end
        LinearAlgebra.:\(s::$XXXMethod, v) = ldiv!(copy(v),s,v)
        
        function $XXXMethod(csr::SparseMatrixCSR{Bi,Tv,Ti},param=nothing) where {Bi,Tv,Ti}
            if csr.m!=csr.n
                error("Matrix must be square")
            end
            myoffset=1-getoffset(Bi)
            csr.rowptr.-=myoffset
            csr.colval.-=myoffset
            s=$XXXMethod(csr.m, csr.rowptr,csr.colval,csr.nzval,tojson(param))
            csr.rowptr.+=myoffset
            csr.colval.+=myoffset
            return s
        end
        
        function $XXXMethod(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},param="") where {Tv,Ti}
            if csc.m!=csc.n
                error("Matrix must be square")
            end
            $XXXMethod(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),param)
        end
    end
end

for XXXMethod in unique([Symbol(instance["XXXMethod"]) for instance in block_instances])
    @eval begin
        LinearAlgebra.ldiv!(s::$XXXMethod, v) = v.=ldiv!(copy(v),s,v)
        function LinearAlgebra.ldiv!(u, s::$XXXMethod, v)
            apply!(s,u,v)
            u
        end
        LinearAlgebra.:\(s::$XXXMethod, v) = ldiv!(copy(v),s,v)
        
        function $XXXMethod(csr::SparseMatrixCSR{Bi,Tv,Ti},blocksize,param=nothing) where {Bi,Tv,Ti}
            if csr.m!=csr.n
                error("Matrix must be square")
            end
            myoffset=1-getoffset(Bi)
            csr.rowptr.-=myoffset
            csr.colval.-=myoffset
            s=$XXXMethod(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,tojson(param))
            csr.rowptr.+=myoffset
            csr.colval.+=myoffset
            return s
        end
        
        function $XXXMethod(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},blocksize,param="") where {Tv,Ti}
            if csc.m!=csc.n
                error("Matrix must be square")
            end
            $XXXMethod(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))),blocksize,param)
        end
    end
end

export AMGSolver,RLXSolver,AMGPrecon,RLXPrecon
export BlockAMGSolver, BlockRLXSolver, BlockAMGPrecon, BlockRLXPrecon

end

