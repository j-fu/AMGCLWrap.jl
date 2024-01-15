# Parameters
```@docs
AMGCLWrap.AbstractAMGCLParams
```

## Iterative solver parameters
See the [Iterative Solvers](https://amgcl.readthedocs.io/en/latest/components/iter_solvers.html#) section
of the [AMGCL documentation](https://amgcl.readthedocs.io/en/latest/index.html).

Instead of one of the solvers below, a named tuple like `(type= "cg", tol=1.0e-10)` can be passed
This also allows to use methods not accessible via one of the structures defined below. 

```@docs
AMGCLWrap.AbstractSolver
AMGCLWrap.BICGStabSolver
AMGCLWrap.CGSolver
AMGCLWrap.GMRESSolver
```

## Relaxation parameters
See the [Relaxation](https://amgcl.readthedocs.io/en/latest/components/relaxation.html) section
of the [AMGCL documentation](https://amgcl.readthedocs.io/en/latest/index.html).

Instead of one of the strategies below, a named tuple like `(type= "damped_jacobi", damping=0.72)` can be passed
This also allows to use methods not accessible via one of the structures defined below. 


```@docs
AMGCLWrap.AbstractRelaxation
AMGCLWrap.ILU0Relaxation
AMGCLWrap.SPAI0Relaxation
```

## Coarsening strategies
See the [Coarsening Strategies](https://amgcl.readthedocs.io/en/latest/components/coarsening.html) section
of the [AMGCL documentation](https://amgcl.readthedocs.io/en/latest/index.html).


Instead of one of the strategies below, a named tuple like `(type= "smoothed_aggregation", relax=1.0)` can be passed
This also allows to use methods not accessible via one of the structures defined below. 

```@docs
AMGCLWrap.AbstractCoarsening
AMGCLWrap.RugeStubenCoarsening
AMGCLWrap.SmoothedAggregationCoarsening
```
