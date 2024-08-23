using AMGCLWrap
using Test, LinearAlgebra, SparseArrays
using Krylov
using LinearSolve
using ILUZero
using Random
import ExplicitImports, Aqua

@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(AMGCLWrap) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(AMGCLWrap) === nothing
end

@testset "Aqua" begin
    Aqua.test_all(AMGCLWrap)
end

if isdefined(Docs,:undocumented_names) # >=1.11
    @testset "UndocumentedNames" begin
        @test isempty(Docs.undocumented_names(AMGCLWrap))
    end
end



A ⊕ B = kron(I(size(B, 1)), A) + kron(B, I(size(A, 1)))

function lattice(n; Tv = Float64)
    d = fill(2 * one(Tv), n)
    d[1] = one(Tv)
    d[end] = one(Tv)
    spdiagm(1 => -ones(Tv, n - 1), 0 => d, -1 => -ones(Tv, n - 1))
end

lattice(L...; Tv = Float64) = lattice(L[1]; Tv) ⊕ lattice(L[2:end]...; Tv)

function dlattice(dim, N; Tv = Float64, Ti = Int64, dd = 1.0e-2, skew=0.0)
    n = N^(1 / dim) |> ceil |> Int
    D=Diagonal(rand(1-skew:0.001:1.0+skew,n^dim))
    s=SparseMatrixCSC{Tv, Ti}(lattice([n for i = 1:dim]...; Tv)*D + Tv(dd) * I)
end;

function iterate(A, f, M)
    u, stats = Krylov.cg(A, f; M, ldiv = true, rtol = 1.0e-12)
    u
end

function test_amg(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    f = A * u0
    amg = AMGSolver(A; blocksize = bsize)
    u = amg \ f
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
end

function test_rlx(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    f = A * u0
    rlx = RLXSolver(A;
                    blocksize = bsize,
                    param = (solver = (tol = 1.0e-12, type = "bicgstab"), precond = (type = "ilu0",)),)
    u = rlx \ f
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
end

function test_amgprecon(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    f = A * u0
    amg = AMGPrecon(A; blocksize = bsize)
    u = iterate(A, f, amg)
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
    true
end



function test_rlxprecon(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    f = A * u0
    rlx = RLXPrecon(A; blocksize = bsize, param = (type = "damped_jacobi",))
    u = iterate(A, f, rlx)
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
    true
end

function test_err(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    f = A * u0
    rlx = RLXSolver(A;
                    blocksize = bsize,
                    param = (solver = (tol = 1.0e-12, type = "bicgstab"), precond = (type = "ilu0x",)),)
    if error_state(rlx) != 0
        println("error catched")
        return true
    end
    return false
end







const NTest = 10000

@testset "catch error" begin
    @test_throws ErrorException test_err(Int64, 1, NTest)
end

@testset "blocksizes" begin
    @test blocksize_instantiated(1)
    @test blocksize_instantiated(2)
    @test !blocksize_instantiated(100)
end

@testset "struct API" begin
    A = dlattice(3, NTest)
    u0 = rand(size(A, 1))
    f = A * u0
    @test isa(AMGSolver(A; coarsening = AMGCLWrap.RugeStubenCoarsening()) \ f, Vector)
    @test isa(RLXSolver(A; solver = AMGCLWrap.CGSolver()) \ f, Vector)
    @test isa(RLXSolver(A; solver = (type = "cg",)) \ f, Vector)
end

for Ti in [Int32, Int64]
    @testset "AMGSolver, $Ti" begin
        @test test_amg(Ti, 1, NTest)
        @test test_amg(Ti, 2, NTest)
        @test test_amg(Ti, 3, NTest)
    end

    @testset "RLXSolver, $Ti" begin
        @test test_rlx(Ti, 1, NTest)
        @test test_rlx(Ti, 2, NTest)
        @test test_rlx(Ti, 3, NTest)
    end

    @testset "AMGPrecon, $Ti" begin
        @test test_amgprecon(Ti, 1, NTest)
        @test test_amgprecon(Ti, 2, NTest)
        @test test_amgprecon(Ti, 3, NTest)
    end

    @testset "RLXPrecon, $Ti" begin
        @test test_rlxprecon(Ti, 1, NTest)
        @test test_rlxprecon(Ti, 2, NTest)
        @test test_rlxprecon(Ti, 3, NTest)
    end

    @testset "blocksize 2, $Ti" begin
        @test test_amg(Ti, 3, NTest, 2)
        @test test_rlx(Ti, 3, NTest, 2)
        @test test_amgprecon(Ti, 3, NTest, 2)
        @test test_rlxprecon(Ti, 3, NTest, 2)
    end
end

function tskew(;skew=0.0)
    Random.seed!(4321)
    dd=0.01
    A = dlattice(3, 100000;skew, dd)
    n=size(A, 1)
    u0 = rand(n)
    f = A * u0
    @show sum(A) ≈ n*dd
    rlx = RLXPrecon(A; param = (type = "ilu0",))
    u1, stats1 = Krylov.bicgstab(A, f; M=rlx, ldiv = true, rtol = 1.0e-12)
    iluz=ilu0(A)
    u2, stats2 = Krylov.bicgstab(A, f; M=iluz, ldiv = true, rtol = 1.0e-12)
    @show norm(u1 - u2)
    @show abs(stats1.niter-stats2.niter)
    norm(u1 - u2) < 10000 * sqrt(eps(Float64))  && abs(stats1.niter-stats2.niter)<5
    
end

@testset "skew" begin
    @test tskew(skew=0.2)
    @test tskew(skew=0.4)
    @test tskew(skew=0.6)
end


function test_linsolve_amg(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)
    u=solve(prb,AMGSolverAlgorithm())
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
end

function test_linsolve_rlx(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)
    u=solve(prb,RLXSolverAlgorithm(;
                                   param = (solver = (tol = 1.0e-12, type = "bicgstab"), precond = (type = "ilu0",)),))
    
    @show norm(u0 - u)
    norm(u0 - u) < 10 * sqrt(eps(Float64))
end

function test_linsolve_amgprecon(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)

    amg = AMGPrecon(A; blocksize = bsize)
    u = solve(prb,KrylovJL_CG(), Pl=amg)
    @show norm(u0 - u)
    norm(u0 - u) < 1000 * sqrt(eps(Float64))
end


function test_linsolve_amgprecs(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)
    u = solve(prb,KrylovJL_CG(precs=AMGPreconditioner(blocksize = bsize)))
    @show norm(u0 - u)
    norm(u0 - u) < 1000 * sqrt(eps(Float64))
end


function test_linsolve_rlxprecon(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)

    rlx = RLXPrecon(A; blocksize = bsize, precond=(type="ilu0",))
    u = solve(prb,KrylovJL_CG(), Pl=rlx)
    @show norm(u0 - u)
    norm(u0 - u) < 1.0e4 * sqrt(eps(Float64))
end

function test_linsolve_rlxprecs(Ti, dim, n, bsize = 1)
    A = dlattice(dim, n; Ti)
    u0 = rand(size(A, 1))
    prb=LinearProblem(A,A*u0)

    rlx =
    u = solve(prb,KrylovJL_CG(precs= RLXPreconditioner(blocksize = bsize, precond=(type="ilu0",))))
    @show norm(u0 - u)
    norm(u0 - u) < 1.0e4 * sqrt(eps(Float64))
end



for Ti in [Int32, Int64]
    @testset "LinearSolve, $Ti" begin
        @test test_linsolve_amg(Ti, 2, NTest)
        @test test_linsolve_rlx(Ti, 2, NTest)
        @test test_linsolve_amgprecon(Ti, 2, NTest)
        @test test_linsolve_amgprecs(Ti, 2, NTest)
        @test test_linsolve_rlxprecon(Ti, 2, NTest)
        @test test_linsolve_rlxprecs(Ti, 2, NTest)
    end

end
