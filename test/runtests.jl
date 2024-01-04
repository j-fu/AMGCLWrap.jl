using AMGCLWrap
using Test, LinearAlgebra, SparseArrays
using Krylov

A ⊕ B = kron(I(size(B, 1)), A) + kron(B, I(size(A, 1)))
function lattice(n; Tv = Float64)
    d = fill(2 * one(Tv), n)
    d[1] = one(Tv)
    d[end] = one(Tv)
    spdiagm(1 => -ones(Tv, n - 1), 0 => d, -1 => -ones(Tv, n - 1))
end

lattice(L...; Tv = Float64) = lattice(L[1]; Tv) ⊕ lattice(L[2:end]...; Tv)

function dlattice(Ti,dim, N; Tv = Float64, dd = 1.0e-2)
    n = N^(1 / dim) |> ceil |> Int
    SparseMatrixCSC{Float64,Ti}(lattice([n for i in 1:dim]...; Tv) + Tv(dd) * I)
end;


function test_amg(Ti,dim,n)
    A=dlattice(Ti,dim,n)
    @show size(A,1)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGSolver(A)
    u=amg\f
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_rlx(Ti,dim,n)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXSolver(A, (solver=(tol=1.0e-12,type="bicgstab"), precond=(type="ilu0",)))
    u=rlx\f
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_amgprecon(Ti,dim,n)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGPrecon(A)
    u,stats=bicgstab(A,f;M=amg,ldiv=true, rtol=1.0e-12)
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_rlxprecon(Ti,dim,n)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXPrecon(A)
    u,stats=bicgstab(A,f;M=rlx,ldiv=true, rtol=1.0e-14,atol=1.0e-20)
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_amg(Ti,dim,n,bsize)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    @show size(A,1)
    f=A*u0
    amg=BlockAMGSolver(A,bsize)
    u=amg\f
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_rlx(Ti,dim,n,bsize)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=BlockRLXSolver(A, bsize, (solver=(tol=1.0e-12,type="bicgstab"), precond=(type="ilu0",)))
    u=rlx\f
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_amgprecon(Ti,dim,n,bsize)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=BlockAMGPrecon(A,bsize)
    u,stats=bicgstab(A,f;M=amg,ldiv=true, rtol=1.0e-12)
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_rlxprecon(Ti,dim,n,bsize)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=BlockRLXPrecon(A,bsize)
    u,stats=bicgstab(A,f;M=rlx,ldiv=true, rtol=1.0e-14,atol=1.0e-20)
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end



const NTest=10000

if Sys.WORD_SIZE == 64
    Tis=[Int32, Int64]
else
    Tis=[Int32]
end

for Ti in Tis
    
@testset "AMGSolver, $Ti" begin
  @test test_amg(Ti,1,NTest)
  @test test_amg(Ti,2,NTest)
  @test test_amg(Ti,3,NTest)
end

@testset "RLXSolver, $Ti" begin
  @test test_rlx(Ti,1,NTest)
  @test test_rlx(Ti,2,NTest)
  @test test_rlx(Ti,3,NTest)
end

@testset "AMGPrecon, $Ti" begin
  @test test_amgprecon(Ti,1,NTest)
  @test test_amgprecon(Ti,2,NTest)
  @test test_amgprecon(Ti,3,NTest)
end

@testset "RLXPrecon, $Ti" begin
  @test test_rlxprecon(Ti,1,NTest)
  @test test_rlxprecon(Ti,2,NTest)
  @test test_rlxprecon(Ti,3,NTest)
end
@testset "blocksize 2, $Ti" begin
 @test test_amg(Ti,3,NTest,2)
 @test test_rlx(Ti,3,NTest,2)
 @test test_amgprecon(Ti,3,NTest,2)
 @test test_rlxprecon(Ti,3,NTest,2)
end

end