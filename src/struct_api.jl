abstract type AbstractAMGCLParams end


abstract type AbstractCoarsening <: AbstractAMGCLParams end
const AbstractCoarseningOrNamedTuple = Union{AbstractCoarsening, NamedTuple}

Base.@kwdef struct SmoothedAggregationCoarsening <: AbstractCoarsening
    type::String="smoothed_aggregation"
    relax::Float64=1.0
end

Base.@kwdef struct RugeStubenCoarsening <: AbstractCoarsening
    type::String="ruge_stuben"
    eps_strong::Float64=0.25
    do_trunc::Bool=true
    eps_trunc::Float64=0.2    
end

abstract type AbstractRelaxation <: AbstractAMGCLParams end
const AbstractRelaxationOrNamedTuple = Union{AbstractRelaxation, NamedTuple}

Base.@kwdef struct SPAI0Relaxation <: AbstractRelaxation
   type::String="spai0"       
end

Base.@kwdef struct ILU0Relaxation <: AbstractRelaxation
   type::String="ilu0"       
end

abstract type AbstractSolver <: AbstractAMGCLParams end
const AbstractSolverOrNamedTuple = Union{AbstractSolver, NamedTuple}


Base.@kwdef struct SolverControl
    tol::Float64=1.0e-10
    abstol::Float64=floatmin(Float64)
    maxiter::Int=100
    verbose::Bool=false
end

@composite Base.@kwdef struct BICGStabSolver <: AbstractSolver
    type::String="bicgstab"
    pside::String="right"
    SolverControl...        
end

@composite Base.@kwdef struct GMRESSolver <: AbstractSolver
    type::String="gmres"
    M::Int=30    
    pside::String="right"
    SolverControl...        
end

@composite Base.@kwdef struct CGSolver <: AbstractSolver
    type::String="cg"
    SolverControl...        
end
    

Base.show(io::IO, p::AbstractAMGCLParams)= JSON3.pretty(io,p)

function AMGSolver(A::AbstractMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   coarsening::AbstractCoarseningOrNamedTuple=SmoothedAggregationCoarsening(),
                   relaxation::AbstractRelaxationOrNamedTuple=SPAI0Relaxation(),
                   solver::AbstractSolverOrNamedTuple=BICGStabSolver(;verbose))
    if param==nothing
            param=(solver=solver, precond=(coarsening=coarsening,relax=relaxation))
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGSolver(A,param;blocksize)
end
    
function RLXSolver(A::AbstractMatrix;
                   param=nothing,
                   verbose::Bool=true,
                   blocksize::Int=1,
                   precond::AbstractRelaxationOrNamedTuple=SPAI0Relaxation(),
                   solver::AbstractSolverOrNamedTuple=BICGStabSolver(;verbose))
    if param==nothing
        param=(solver=solver, precond=precond)
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXSolver(A,param;blocksize)
end
    
function AMGPrecon(A::AbstractMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   coarsening::AbstractCoarseningOrNamedTuple=SmoothedAggregationCoarsening(),
                   relax::AbstractRelaxationOrNamedTuple=SPAI0Relaxation())
    if param==nothing
        param=(coarsening=coarsening,relax=relax)
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGPrecon(A,param;blocksize)
end
    
function RLXPrecon(A::AbstractMatrix;
                   param=nothing,
                   verbose::Bool=false,
                   blocksize::Int=1,
                   precond::AbstractRelaxationOrNamedTuple=SPAI0Relaxation())
    if param==nothing
         param=relaxation
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXPrecon(A,param;blocksize)
end
    

    
