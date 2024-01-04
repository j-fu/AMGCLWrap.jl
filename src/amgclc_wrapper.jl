# 
# Iteration info returned by solver application methods
#
struct AMGCLInfo
    iters::Cint
    residual::Cdouble
end


# Create json string from Tuple, Dict, String or Nothing
tojson(param)=JSON3.write(param)
tojson(::Nothing)=""
tojson(s::String)=s

#
# Operators to be created 
#
const operators=[:AMGSolver,:RLXSolver,:AMGPrecon,:RLXPrecon]

const docs=Dict(
    "AMGCLWrap.AMGSolver"=>"Create Algebraic multigrid preconditioned Krylov subspace solver with `ldiv!` and `\\` methods solving the matrix system.",
    "AMGCLWrap.RLXSolver"=>"Create single level relaxation preconditioned Krylov subspace solver with `ldiv!` and `\\` methods solving the matrix system.",
    "AMGCLWrap.AMGPrecon"=>"Create algebraic multigrid preconditioner with `ldiv!` and `\\` methods solving the preconditioning system.",
    "AMGCLWrap.RLXPrecon"=>"Create single level relaxation preconditioner with `ldiv!` and `\\` methods solving the preconditioning system."
)


#
# Dicts to support type mangling 
#
const dtypedict=Dict(
    "RealChar"=> "D",
    "JTv"=>"Float64",
    "CTv"=>"Cdouble")

const ltypedict=Dict(
    "IntChar"=> "L",
    "JTi"=>"Int64",
    "CTi"=>"Clong")

const itypedict=Dict(
    "IntChar"=> "I",
    "JTi"=>"Int32",
    "CTi"=>"Cint")


if Sys.WORD_SIZE == 64
    const idxtypedicts=[ltypedict,itypedict]
else
    const idxtypedict=[itypedict]
end    

#
# Dicts which describe how to wrap operator names in C
#
const operatordicts=[]

for method in operators
    for idxtypedict in idxtypedicts
        for numtypedict in [dtypedict]
            TvTi=numtypedict["RealChar"]*idxtypedict["IntChar"]
            dict=Dict(
                "Operator" => String(method),
                "TvTi" => TvTi,
                "JTv"  => numtypedict["JTv"],
                "JTi"  => idxtypedict["JTi"],
                "CTv"  => numtypedict["CTv"],
                "CTi"  => idxtypedict["CTi"])
            push!(operatordicts,dict)
        end
    end
end


abstract type AbstractAMGCLOperator end
#
# Define operator types
#
for Operator in operators
    @eval begin
        mutable struct $Operator{Tv,Ti} <: AbstractAMGCLOperator
            handle::Ptr{Cvoid}
            blocksize::Cint
        end
    end
end

#
# Define julia methods wrapping C functions
#
for operatordict in operatordicts
    Operator=Symbol(operatordict["Operator"])
    TvTiOperator=operatordict["TvTi"]*operatordict["Operator"]
    amgclcTvTiOperatorCreate="amgclc$(TvTiOperator)Create"
    amgclcTvTiOperatorDestroy="amgclc$(TvTiOperator)Destroy"
    amgclcTvTiOperatorApply="amgclc$(TvTiOperator)Apply"
    JTv=Symbol(operatordict["JTv"])
    JTi=Symbol(operatordict["JTi"])
    CTv=Symbol(operatordict["CTv"])
    CTi=Symbol(operatordict["CTi"])

    @eval begin
        #
        # Constructor from bunch of arrays
        #
        function $Operator(n, ia::Vector{Ti}, ja::Vector{Ti}, a::Vector{Tv}, blocksize,param::String) where {Tv<:$JTv,Ti<:$JTi}
            this=ccall(($amgclcTvTiOperatorCreate,libamgcl_c),
                       $Operator{$JTv,$JTi},
                       (Cint, Ptr{$CTi}, Ptr{$CTi},Ptr{$CTv},Cint,Cstring),
                       n, ia,ja, a, blocksize, param);
            finalizer(s->ccall(($amgclcTvTiOperatorDestroy,libamgcl_c),
                               Cvoid,
                               ($Operator{$JTv, $JTi},),
                               s),this)
            this
        end

        #
        # Solve matrix/preconditioning system
        #
        function apply!(s::$Operator{Tv,Ti}, sol::Vector{Tv}, rhs::Vector{Tv})  where {Tv<:$JTv,Ti<:$JTi}
            i=ccall(($amgclcTvTiOperatorApply,libamgcl_c),
                    AMGCLInfo,
                    ($Operator{$JTv, $JTi},Ptr{$CTv}, Ptr{$CTv}),
                    s,sol,rhs)
            (iters=i.iters, residual=i.residual)
        end
    end
end

#
# Linear Algebra solution methods
#
function LinearAlgebra.ldiv!(u, s::AbstractAMGCLOperator, v)
    apply!(s,u,v)
    u
end
LinearAlgebra.ldiv!(s::AbstractAMGCLOperator, v) = v.=ldiv!(copy(v),s,v)
LinearAlgebra.:\(s::AbstractAMGCLOperator, v) = ldiv!(copy(v),s,v)


#
# Constructors from sparse matrices
#
for Operator in operators
    @eval begin
        @doc """
           $($Operator)(sparsematrix; blocksize=1, param=nothing)

           $(docs[string($Operator)])

           Input: 
            - `sparsematrix`: `SparseArrays.AbstractSparseMatrixCSC` or `SparseMatricesCSR.SparseMatrixCSR`. 
            - `blocksize`: If blocksize >1, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of        `blocksize x blocksize` static matrices. Block sizes 1...8 are instantiated.
            - `params`: Any object (e.g. Tuple, Dict or JSON string) which can be turned into a JSON string by `JSON3.write`. If `params` is an emtpy string or `nothing` a default value is used.
         """
        function $Operator(csr::SparseMatrixCSR{Bi,Tv,Ti}; blocksize=1, param=nothing) where {Bi,Tv,Ti}
            if csr.m!=csr.n
                error("Matrix must be square")
            end
            myoffset=1-getoffset(Bi)
            csr.rowptr.-=myoffset
            csr.colval.-=myoffset
            s=$Operator(csr.m, csr.rowptr,csr.colval,csr.nzval,blocksize,tojson(param))
            csr.rowptr.+=myoffset
            csr.colval.+=myoffset
            return s
        end
        function $Operator(csc::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}; blocksize=1,param=nothing) where {Tv,Ti}
            $Operator(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc)));blocksize,param)
        end
    end
end

