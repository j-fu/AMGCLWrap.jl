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

function dlattice(dim, N; Tv = Float64, dd = 1.0e-2)
    n = N^(1 / dim) |> ceil |> Int
    lattice([n for i in 1:dim]...; Tv) + Tv(dd) * I
end;


function test_amg(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGSolver(A)
    u=zeros(size(A,1))
    info=AMGCLWrap.apply!(amg,u,f)
    @show info
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_rlx(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXSolver(A, (solver=(tol=1.0e-12,type="bicgstab"), precond=(type="ilu0",)))
    u=zeros(size(A,1))
    info=AMGCLWrap.apply!(rlx,u,f)
    @show info
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_amg2(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGSolver(A)
    u=amg\f
    norm(u0-u)<sqrt(eps(Float64))
end

function test_rlx2(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXSolver(A, (solver=(tol=1.0e-12,type="bicgstab"), precond=(type="ilu0",)))
    u=rlx\f
    norm(u0-u)<sqrt(eps(Float64))
end

function test_amgprecon(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGPrecon(A)
    u,stats=bicgstab(A,f;M=amg,ldiv=true, rtol=1.0e-12)
    @show stats
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end

function test_rlxprecon(dim,n)
    A=dlattice(dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXPrecon(A)
    u,stats=bicgstab(A,f;M=rlx,ldiv=true, rtol=1.0e-14,atol=1.0e-20)
    @show stats
    @show norm(u0-u)
    norm(u0-u)<sqrt(eps(Float64))
end


@testset "AMGSolver" begin

  @test test_amg(1,10000)
  @test test_amg(2,10000)
  @test test_amg(3,10000)
    
  @test test_amg2(1,10000)
  @test test_amg2(2,10000)
  @test test_amg2(3,10000)
    
end

@testset "RLXSolver" begin

  @test test_rlx(1,10000)
  @test test_rlx(2,10000)
  @test test_rlx(3,10000)
    
  @test test_rlx2(1,10000)
  @test test_rlx2(2,10000)
  @test test_rlx2(3,10000)
end

@testset "AMGPrecon" begin

  @test test_amgprecon(1,10000)
  @test test_amgprecon(2,10000)
  @test test_amgprecon(3,10000)
    
    
end

@testset "RLXPrecon" begin

  @test test_rlxprecon(1,10000)
  @test test_rlxprecon(2,10000)
  @test test_rlxprecon(3,10000)
    
    
end
