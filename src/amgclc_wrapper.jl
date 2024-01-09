
# iteration/solver info returned by amgcl_c
struct AMGCLInfo
    iters::Cint
    residual::Cdouble
    error_state::Cint
end

# Create json string from Tuple, Dict, String or Nothing
tojson(param) = JSON3.write(param)
tojson(::Nothing) = ""
tojson(s::String) = s

"""
    blocksize_instantiated(blocksize)::Bool

Check if given blocksize has been instantiated.
"""
blocksize_instantiated(blocksize) = ccall((:amgclcBlocksizeInstantiated, libamgcl_c), Cint, (Cint,), blocksize) == 1

#
# Operators to be created 
#
const operators = [:AMGSolver, :RLXSolver, :AMGPrecon, :RLXPrecon]

const docs = Dict("AMGCLWrap.AMGSolver" => "Create Algebraic multigrid preconditioned Krylov subspace solver with `ldiv!` and `\\` methods solving the matrix system.",
                  "AMGCLWrap.RLXSolver" => "Create single level relaxation preconditioned Krylov subspace solver with `ldiv!` and `\\` methods solving the matrix system.",
                  "AMGCLWrap.AMGPrecon" => "Create algebraic multigrid preconditioner with `ldiv!` and `\\` methods solving the preconditioning system.",
                  "AMGCLWrap.RLXPrecon" => "Create single level relaxation preconditioner with `ldiv!` and `\\` methods solving the preconditioning system.")

#
# Dicts to support type mangling 
#
const dtypedict = Dict("RealChar" => "D",
                       "JTv" => "Float64",
                       "CTv" => "Cdouble")

#
# amgcl_c guarantees that "DI" maps to Int32 and "DL" to Int64
#
const ltypedict = Dict("IntChar" => "L",
                       "JTi" => "Int64",
                       "CTi" => "Int64")

const itypedict = Dict("IntChar" => "I",
                       "JTi" => "Int32",
                       "CTi" => "Int32")

const idxtypedicts = [ltypedict, itypedict]

const numtypedicts = [dtypedict]
#
# Dicts which describe how to wrap operator names in C
#
const operatordicts = []

for operator in operators
    for idxtypedict in idxtypedicts
        for numtypedict in numtypedicts
            TvTi = numtypedict["RealChar"] * idxtypedict["IntChar"]
            dict = Dict("Operator" => String(operator),
                        "TvTi" => TvTi,
                        "JTv" => numtypedict["JTv"],
                        "JTi" => idxtypedict["JTi"],
                        "CTv" => numtypedict["CTv"],
                        "CTi" => idxtypedict["CTi"])
            push!(operatordicts, dict)
        end
    end
end

abstract type AbstractAMGCLOperator end

"""
    error_state(operator)::Int
Return error state of operator. Ok if 0.
"""
error_state(o::AbstractAMGCLOperator) = o.error_state

#
# Define operator types
#
for Operator in operators
    @eval begin
        mutable struct $Operator{Tv, Ti} <: AbstractAMGCLOperator
            handle::Ptr{Cvoid}
            blocksize::Cint
            error_state::Cint
        end
    end
end

issolver(::Any) = false
issolver(::AMGSolver) = true
issolver(::RLXSolver) = true

#
# Define julia methods wrapping C functions
#
for operatordict in operatordicts
    Operator = Symbol(operatordict["Operator"])
    TvTiOperator = operatordict["TvTi"] * operatordict["Operator"]
    amgclcTvTiOperatorCreate = "amgclc$(TvTiOperator)Create"
    amgclcTvTiOperatorDestroy = "amgclc$(TvTiOperator)Destroy"
    amgclcTvTiOperatorApply = "amgclc$(TvTiOperator)Apply"
    JTv = Symbol(operatordict["JTv"])
    JTi = Symbol(operatordict["JTi"])
    CTv = Symbol(operatordict["CTv"])
    CTi = Symbol(operatordict["CTi"])

    @eval begin
        function finalize!(operator::$Operator{$JTv, $JTi})
            ccall(($amgclcTvTiOperatorDestroy, libamgcl_c),
                  Cvoid,
                  ($Operator{$JTv, $JTi},),
                  operator)
        end

        #
        # Constructor from bunch of arrays
        #
        function $Operator(n, ia::Vector{$JTi}, ja::Vector{$JTi}, a::Vector{$JTv}, blocksize, param::String)
            this = ccall(($amgclcTvTiOperatorCreate, libamgcl_c),
                         $Operator{$JTv, $JTi},
                         (Cint, Ptr{$CTi}, Ptr{$CTi}, Ptr{$CTv}, Cint, Cstring),
                         n, ia, ja, a, blocksize, param)
            finalizer(finalize!, this)
            this
        end

        function apply!(operator::$Operator{$JTv, $JTi}, sol::Vector{$JTv}, rhs::Vector{$JTv})
            info = ccall(($amgclcTvTiOperatorApply, libamgcl_c),
                         AMGCLInfo,
                         ($Operator{$JTv, $JTi}, Ptr{$CTv}, Ptr{$CTv}),
                         operator, sol, rhs)
            return (iters = info.iters, residual = info.residual, error_state = info.error_state)
        end
    end
end

#
# Linear Algebra solution methods
#
function LinearAlgebra.ldiv!(u, operator::AbstractAMGCLOperator, v)
    apply!(operator, u, v)
    u
end
function LinearAlgebra.ldiv!(operator::AbstractAMGCLOperator, v)
    u = ldiv!(copy(v), operator, v)
    v .= u
    nothing
end

LinearAlgebra.:\(operator::AbstractAMGCLOperator, v) = ldiv!(copy(v), operator, v)

#
# Constructors from sparse matrices
#
for Operator in operators
    @eval begin
        @doc """

             $($Operator)(sparsematrix::AbstractSparseMatrix; 
                          blocksize=1, 
                          param=nothing)


           $(docs[string($Operator)])

           Parameters:
            - `sparsematrix`: `SparseArrays.AbstractSparseMatrixCSC` or `SparseMatricesCSR.SparseMatrixCSR`. 
            - `blocksize`: If blocksize >1, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of        `blocksize x blocksize` static matrices. Block sizes 1...8 are instantiated.
            - `param`: Any object (e.g. Tuple, Dict or JSON string) which can be turned into a JSON string by `JSON3.write`. If `params` is an emtpy string or `nothing` a default value is used.
         """
        function $Operator(csr::SparseMatrixCSR{Bi, Tv, Ti}, param; blocksize = 1) where {Bi, Tv, Ti}
            if csr.m != csr.n
                error("Matrix must be square")
            end
            myoffset = 1 - getoffset(Bi)
            csr.rowptr .-= myoffset
            csr.colval .-= myoffset
            operator = $Operator(csr.m, csr.rowptr, csr.colval, csr.nzval, blocksize, tojson(param))
            csr.rowptr .+= myoffset
            csr.colval .+= myoffset
            return operator
        end

        function $Operator(csc::SparseArrays.AbstractSparseMatrixCSC{Tv, Ti}, param; blocksize = 1) where {Tv, Ti}
            $Operator(SparseMatrixCSR{1}(transpose(SparseMatrixCSC(csc))), param; blocksize)
        end
    end
end
