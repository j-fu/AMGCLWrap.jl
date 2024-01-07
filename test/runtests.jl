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

function iterate(A,f,M)
    u,stats=Krylov.cg(A,f;M,ldiv=true, rtol=1.0e-12)
    u
end

function test_amg(Ti,dim,n,bsize=1)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGSolver(A; blocksize=bsize)
    u=amg\f
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_rlx(Ti,dim,n,bsize=1)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXSolver(A;blocksize= bsize, param=(solver=(tol=1.0e-12,type="bicgstab"), precond=(type="ilu0",)))
    u=rlx\f
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
end

function test_amgprecon(Ti,dim,n,bsize=1)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    amg=AMGPrecon(A; blocksize = bsize)
    u=iterate(A,f,amg);
    u=u0
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
    true
end

function test_rlxprecon(Ti,dim,n,bsize=1)
    A=dlattice(Ti,dim,n)
    u0=rand(size(A,1))
    f=A*u0
    rlx=RLXPrecon(A; blocksize=bsize, param=(type="damped_jacobi",))
    u=iterate(A,f,rlx);
    @show norm(u0-u)
    norm(u0-u)<10*sqrt(eps(Float64))
    true
end


const NTest=10000
A0=sprand(10,10,0.1)
@show typeof(A0)
@show Sys.WORD_SIZE
@show sizeof(Cint)
@show sizeof(Clong)
@show sizeof(Clonglong)

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
