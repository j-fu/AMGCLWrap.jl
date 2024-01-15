################################################################################################
# Docstring snippets

const matrixparam="""
                  - `sparsematrix`: `SparseArrays.AbstractSparseMatrixCSC` or `SparseMatricesCSR.SparseMatrixCSR`. 
                  """

const stdparams = """
                  - `sparsematrix`: `SparseArrays.AbstractSparseMatrixCSC` or `SparseMatricesCSR.SparseMatrixCSR`. 
                  - `blocksize`: If blocksize >1, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of `blocksize x blocksize` static matrices. Block sizes 1...8 are instantiated.
                  - `verbose`: if true, print generated JSON string passed to amgcl.
                  - `param:`   Ignored if `nothing` (default). Otherwise, any object (e.g. Tuple, Dict or JSON string) which can be turned into a JSON string by `JSON3.write`.
                  """
const amgsolverparams= """
                       - `coarsening`: One of the  [Coarsening strategies](@ref)
                       - `relax`: One of the [Relaxation strategies](@ref)
                       - `solver`: One of the [Iterative solver strategies](@ref)
                       """    
const rlxsolverparams= """
                       - `precond`: One of the [Relaxation strategies](@ref) seen as preconditioner
                       - `solver`: One of the [Iterative solver strategies](@ref)
                       """    

#################################################################################################
# AMG Solver

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
$(matrixparam)
$(stdparams)
$(amgsolverparams)
"""
function AMGSolver(sparsematrix::AbstractSparseMatrix;
                   param = nothing,
                   verbose::Bool = false,
                   blocksize::Int = 1,
                   coarsening::Union{AbstractCoarsening, NamedTuple} = SmoothedAggregationCoarsening(),
                   relax::Union{AbstractRelaxation, NamedTuple} = SPAI0Relaxation(),
                   solver::Union{AbstractSolver, NamedTuple} = BICGStabSolver(; verbose))
    if param == nothing
        param = (solver = solver, precond = (coarsening = coarsening, relax = relax))
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGSolver(sparsematrix, param; blocksize)
end

@static if VERSION >=v"1.9.0"

Base.@kwdef mutable struct AMGSolverAlgorithmData
    param = nothing
    verbose::Bool = false
    blocksize::Int = 1
    coarsening::Union{AbstractCoarsening, NamedTuple} = SmoothedAggregationCoarsening()
    relax::Union{AbstractRelaxation, NamedTuple} = SPAI0Relaxation()
    solver::Union{AbstractSolver, NamedTuple} = BICGStabSolver(; verbose)
    instance=nothing
end

function (data::AMGSolverAlgorithmData)(A, b, u, p, newA, Pl, Pr, solverdata; verbose = true, kwargs...)
    (;param,verbose,blocksize,solver,coarsening,relax)=data
    if data.instance==nothing || newA
        data.instance=AMGSolver(A;param,verbose,blocksize,solver,relax,coarsening)
    end
    ldiv!(u,data.instance,b)
end

    
end # @static if

"""
    AMGSolverAlgorithm(;blocksize::Int=1,
                        param=nothing,
                        verbose::Bool=false,
                        coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
                        relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation(),
                        solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))

Algebraic multigrid preconditioned Krylov subspace solver algorithm for LinearSolve.jk
Parameters:
$(stdparams)
$(amgsolverparams)

!!! compat
    Only available for Julia version >=1.9
"""
function AMGSolverAlgorithm end

#################################################################################################
# Relaxation  Solver
"""
    RLXSolver(sparsematrix::AbstractSparseMatrix;
              blocksize::Int=1,
              param=nothing,
              verbose::Bool=false,
              precond::Union{AbstractRelaxation, NamedTuple}=ILU0Relaxation(),
              solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))
$(docs["AMGCLWrap.RLXSolver"])

Parameters:
$(matrixparam)
$(stdparams)
$(rlxsolverparams)
"""
function RLXSolver(sparsematrix::AbstractSparseMatrix;
                   param = nothing,
                   verbose::Bool = false,
                   blocksize::Int = 1,
                   precond::Union{AbstractRelaxation, NamedTuple} = ILU0Relaxation(),
                   solver::Union{AbstractSolver, NamedTuple} = BICGStabSolver())
    if param == nothing
        param = (solver = solver, precond = precond)
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXSolver(sparsematrix, param; blocksize)
end


@static if VERSION >=v"1.9.0"

Base.@kwdef mutable struct RLXSolverAlgorithmData
    param = nothing
    verbose::Bool = false
    blocksize::Int = 1
    precond::Union{AbstractRelaxation, NamedTuple} = SPAI0Relaxation()
    solver::Union{AbstractSolver, NamedTuple} = BICGStabSolver(; verbose)
    instance=nothing
end

function (data::RLXSolverAlgorithmData)(A, b, u, p, newA, Pl, Pr, solverdata; verbose = true, kwargs...)
    (;param,verbose,blocksize,solver,precond)=data
    if data.instance==nothing || newA
        data.instance=RLXSolver(A;param,verbose,blocksize,solver,precond)
    end
    ldiv!(u,data.instance,b)
end

end # @static if



"""
    RLXSolverAlgorithm(;blocksize::Int=1,
                        param=nothing,
                        verbose::Bool=false,
                        precond::Union{AbstractRelaxation, NamedTuple}=ILU0Relaxation(),
                        solver::Union{AbstractSolver, NamedTuple}=BICGStabSolver(;verbose))

Algebraic multigrid preconditioned Krylov subspace solver algorithm for LinearSolve.jk
Parameters:
$(stdparams)
$(rlxsolverparams)

!!! compat
    Only available for Julia version >=1.9
"""
function RLXSolverAlgorithm end

    
#################################################################################################
# AMG Preconditioner

"""
   AMGPrecon(sparsematrix::AbstractSparseMatrix;
             blocksize::Int=1,
             param=nothing,
             verbose::Bool=false,
             coarsening::Union{AbstractCoarsening, NamedTuple}=SmoothedAggregationCoarsening(),
             relax::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
$(docs["AMGCLWrap.AMGPrecon"])

Parameters:
$(matrixparam)
$(stdparams)
- `coarsening`: A [coarsening strategy](#Coarsening-strategies)
- `relax`: A [relaxation method](#Relaxation/Preconditioner-parameters)
"""
function AMGPrecon(sparsematrix::AbstractSparseMatrix;
                   param = nothing,
                   verbose::Bool = false,
                   blocksize::Int = 1,
                   coarsening::Union{AbstractCoarsening, NamedTuple} = SmoothedAggregationCoarsening(),
                   relax::Union{AbstractRelaxation, NamedTuple} = SPAI0Relaxation())
    if param == nothing
        param = (coarsening = coarsening, relax = relax)
    end
    if verbose
        JSON3.pretty(param)
    end
    AMGPrecon(sparsematrix, param; blocksize)
end

"""
    RLXPrecon(sparsematrix::AbstractSparseMatrix;
              blocksize::Int=1,
              param=nothing,
              verbose::Bool=false,
              precond::Union{AbstractRelaxation, NamedTuple}=SPAI0Relaxation())
$(docs["AMGCLWrap.RLXPrecon"])

Parameters:
$(matrixparam)
$(stdparams)
- `precond`: A [preconditioned method](#Relaxation/Preconditioner-parameters)
"""

#################################################################################################
# Relaxation Preconditioner

function RLXPrecon(sparsematrix::AbstractSparseMatrix;
                   param = nothing,
                   verbose::Bool = false,
                   blocksize::Int = 1,
                   precond::Union{AbstractRelaxation, NamedTuple} = SPAI0Relaxation())
    if param == nothing
        param = precond
    end
    if verbose
        JSON3.pretty(param)
    end
    RLXPrecon(sparsematrix, param; blocksize)
end

