module AMGCLWrapLinearSolveExt
import AMGCLWrap: AMGSolverAlgorithmData, AMGSolverAlgorithm
import AMGCLWrap: RLXSolverAlgorithmData, RLXSolverAlgorithm
using LinearSolve: LinearSolveFunction

AMGSolverAlgorithm(;kwargs...)=LinearSolveFunction(AMGSolverAlgorithmData(;kwargs...))
RLXSolverAlgorithm(;kwargs...)=LinearSolveFunction(RLXSolverAlgorithmData(;kwargs...))


end
