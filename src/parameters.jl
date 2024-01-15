#####################################################################
# Apex parameter type

"""
    $(TYPEDEF)

Abstract parameter type.
"""
abstract type AbstractAMGCLParams end

Base.show(io::IO, p::AbstractAMGCLParams) = JSON3.pretty(io, p)

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
    tol::Float64 = 1.0e-10
    abstol::Float64 = floatmin(Float64)
    maxiter::Int = 100
    verbose::Bool = false
end

"""
    $(TYPEDEF)

BICGStab solver
$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct BICGStabSolver <: AbstractSolver
    type::String = "bicgstab"
    pside::String = "right"
    SolverControl...
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct GMRESSolver <: AbstractSolver
    type::String = "gmres"
    M::Int = 30
    pside::String = "right"
    SolverControl...
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
@composite Base.@kwdef struct CGSolver <: AbstractSolver
    type::String = "cg"
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
    type::String = "spai0"
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct ILU0Relaxation <: AbstractRelaxation
    type::String = "ilu0"
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
    type::String = "smoothed_aggregation"
    relax::Float64 = 1.0
end

"""
    $(TYPEDEF)

$(TYPEDFIELDS)
"""
Base.@kwdef struct RugeStubenCoarsening <: AbstractCoarsening
    type::String = "ruge_stuben"
    eps_strong::Float64 = 0.25
    do_trunc::Bool = true
    eps_trunc::Float64 = 0.2
end

