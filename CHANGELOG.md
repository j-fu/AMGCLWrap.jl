# Changelog

## v2.0, 2024-10-23
- Rename Preconditioner to PreconBuilder, closer to the logic of the LinearSolve precs API

## v1.0, 2024-08-23
- Support `precs` API of LinearSolve
- Not really breaking, but transition to genuine semantig versioning


## v0.4, 2024-07-15
- Fix SparseMatrixCSR from CSC for nonsymmetric matrices
- Add Aqua, ExplicitImports tests, ci on apple silicon
- Require julia 1.9

## v0.3, 2024-01-15
- LinearSolve extension, solver support  via a LinearSolveFunction
- Change default for RLXPrecon to ILU0

## v0.2, 2024-01-09
- Finalization of API, exception handling

## v0.1, 2024-01-08
- Initial registered version
