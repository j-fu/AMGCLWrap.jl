#####################################################################
# Apex parameter type

"""
    $(TYPEDEF)

Abstract parameter type.
"""
abstract type AbstractAMGCLParams end

Base.show(io::IO, p::AbstractAMGCLParams)= JSON3.pretty(io,p)

#####################################################################
# Iterative solvers

"""
    $(TYPEDEF)

Abstract solver parameter type.
"""
abstract type AbstractSolver <: AbstractAMGCLParams end

Base.@kwdef struct SolverControl
    """
    Relative tolerance
    """
    tol::Float64=1.0e-10
    abstol::Float64=floatmin(Float64)
    maxiter::Int=100
    verbose::Bool=false
end

"""
    $(TYPEDEF)

BICGStab solver
$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct BICGStabSolver <: AbstractSolver
    type::String="bicgstab"
    pside::String="right"
    SolverControl...        
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct GMRESSolver <: AbstractSolver
    type::String="gmres"
    M::Int=30    
    pside::String="right"
    SolverControl...        
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct CGSolver <: AbstractSolver
    type::String="cg"
    SolverControl...        
end
    
########################################################################
# Relaxation/Preconditioning    

"""
    $(TYPEDEF)

Abstract relaxation parameter type.
"""
abstract type AbstractRelaxation <: AbstractAMGCLParams end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct SPAI0Relaxation <: AbstractRelaxation
   type::String="spai0"       
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct ILU0Relaxation <: AbstractRelaxation
   type::String="ilu0"       
end


#####################################################################
# Coarsening

"""
    $(TYPEDEF)

Abstract coarsening parameter type.
"""
abstract type AbstractCoarsening <: AbstractAMGCLParams end


"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct SmoothedAggregationCoarsening <: AbstractCoarsening
    type::String="smoothed_aggregation"
    relax::Float64=1.0
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct RugeStubenCoarsening <: AbstractCoarsening
    type::String="ruge_stuben"
    eps_strong::Float64=0.25
    do_trunc::Bool=true
    eps_trunc::Float64=0.2    
end
    


const stdparams=
"""
- `sparsematrix`: `SparseArrays.AbstractSparseMatrixCSC` or `SparseMatricesCSR.SparseMatrixCSR`. 
- `blocksize`: If blocksize >1, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of        `blocksize x blocksize` static matrices. Block sizes 1...8 are instantiated.
- `param:`   Ignored if `nothing` (default). Otherwise, any object (e.g. Tuple, Dict or JSON string) which can be turned into a JSON string by `JSON3.write`.
- `verbose`: if true, print generated JSON string passed to amgcl.
"""    

"""
    AMGSolver(sparsematrix::AbstractSparseMatrix;
              blocksize::Int=1,
              param=nothing,
              verbose::Bool=false,
              coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
              relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation(),
              solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))
$(docs["AMGCLWrap.AMGSolver"])

Parameters:
$(stdparams)
- `coarsening`: A [coarsening strategy](#Coarsening-strategies)
- `relax`: A [relaxation method](#Relaxation/Preconditioner-parameters)
- `solver`: An [iterative solver](#Solver-parameters)
"""
function AMGSolver(sparsematrix::AbstractSparseMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
                   relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation(),
                   solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))
    if param==nothing
            param=(solver=solver, precond=(coarsening=coarsening,relax=relax))
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGSolver(sparsematrix,param;blocksize)
end

"""
    RLXSolver(sparsematrix::AbstractSparseMatrix;
              blocksize::Int=1,
              param=nothing,
              verbose::Bool=false,
              precond::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation(),
              solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))
$(docs["AMGCLWrap.RLXSolver"])

Parameters:
$(stdparams)
- `precond`: A [preconditioned method](#Relaxation/Preconditioner-parameters)
- `solver`: An [iterative solver](#Solver-parameters)
"""
function RLXSolver(sparsematrix::AbstractSparseMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   precond::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation(),
                   solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver())
    if param==nothing
        param=(solver=solver, precond=precond)
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXSolver(sparsematrix,param;blocksize)
end
    
"""
   AMGPrecon(sparsematrix::AbstractSparseMatrix;
             blocksize::Int=1,
             param=nothing,
             verbose::Bool=false,
             coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
             relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
$(docs["AMGCLWrap.AMGPrecon"])

Parameters:
$(stdparams)
- `coarsening`: A [coarsening strategy](#Coarsening-strategies)
- `relax`: A [relaxation method](#Relaxation/Preconditioner-parameters)
"""
function AMGPrecon(sparsematrix::AbstractSparseMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
                   relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
    if param==nothing
        param=(coarsening=coarsening,relax=relax)
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGPrecon(sparsematrix,param;blocksize)
end

"""
    RLXPrecon(sparsematrix::AbstractSparseMatrix;
              blocksize::Int=1,
              param=nothing,
              verbose::Bool=false,
              precond::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
$(docs["AMGCLWrap.RLXPrecon"])

Parameters:
$(stdparams)
- `precond`: A [preconditioned method](#Relaxation/Preconditioner-parameters)
"""
function RLXPrecon(sparsematrix::AbstractSparseMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   precond::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
    if param==nothing
         param=relaxation
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXPrecon(sparsematrix,param;blocksize)
end
    

    
