# Some Integration Routines In C #

Just a variety of straight-forward numerica integration routines for real-valued functions on a finite domain.  Single-threaded.


## To include in a project

In your source file include `integrate.h` and call the relevant integration routines.  To build, include `integrate.c` and `interval.c` in your source files. Both require `interval.h` and `integrate.h` to be findable.


## Fixed stencil routines

There are two routines with fixed stencils:

* `trap`: Trapezoid rule, second order convergence.
* `simp`: Simpsons rule, fourth order convergence.

Both have the same call signature:

```c
double trap(double (*f)(double, void *), double xa, double xb, int N,
            void *args, int (*errf)(void *))
```

The arguments are:

* `f`: A pointer to the function to integrate, called as `y = f(x, args)`.
* `xa`, `xb`: Left and right integration bounds
* `N`: Number of evaluations.
* `args`: pointer to the arguments passed to `f`
* `errf`: pointer to the error function. Is called after every evaluation of `f`.  If return value is non-zero, integration is immediately halted.

## Romberg Routine

A simple adaptive Romberg integrator. Exits after error tolerance met, maximum number of evalutions met, or 20 levels have been performed.  Each level internally produces an integal estimate `I` with error `eps`. Integration is terminated when `|eps| < atol + rtol * |I|`.

```c
double romb(double (*f)(double, void *), double xa, double xb, int N,
            double atol, double rtol, void *args, int *Neval, double *eps,
            int verbose, int (*errf)(void *), double *pfa, double *pfb)
```

The arguments are:

* `f`: A pointer to the function to integrate, called as `y = f(x, args)`.
* `xa`, `xb`: Left and right integration bounds
* `N`: Maximum number of evaluations, ignored if < 1. 
* `atol`: Absolute error tolerance
* `rtol`: Relative error tolerance
* `args`: pointer to the arguments passed to `f`
* `errf`: pointer to the error function. Is called after every evaluation
    of `f`.  If return value is non-zero, integration is immediately halted
* `Neval`: pointer to integer. if not NULL, set to the number of evaluations performed
* `eps`: pointer to double. if not NULL, set to the error estimate
* `verbose`: if non-zero, print internal information during run
* `errf`: pointer to the error function. Is called after every evaluation of `f`.  If return value is non-zero, integration is immediately halted.
* `pfa`, `pfb`: pointers to double. If not NULL, used for `f(xa)` and `f(xb)`.

## Adaptive Routines

There are several adaptive routines. Each use a standard adaptive quadrature routine which divides the integration domain into Intervals and calculating the integral in each Interval with a simple stencil. The total integral and error are the sums (in absolute value in case of the error) of contributions from each Interval. At each step of integration, the Interval with the worst error is processed by splitting into two and re-computing the stencil.  Integration terminates when a maximum number of evaluations have been performed, the estimated error is within the tolerance, or an error occurs.

Intervals are stored and sorted in a custom heap `struct Mesh`.  A few versions of both the `struct Interval` and `struct Mesh` are included which store extra data in each interval, to aid in re-using computations.

The adaptive integration routines are:

* `trap_adapt`: Adaptive trapezoid rule.  Slow but simple and robust.
* `simp_adapt`: Adaptive Simpson's rule.
* `gk49_adapt`: Adaptive Gauss-Kronrod: G4 K9.
* `gk715_adapt`: Adaptive Gauss-Kronrod: G7 K15.
* `gk1021_adapt`: Adaptive Gauss-Kronrod: G10 K21.
* `cadre_adapt`: A custom routine based on the CADRE integrator. Begins with a multi-point trapezoid rule in each interval to get a non-linear error estimate and an estimate of the convergence rate. While the convergence rate is not second order (implying we are under-resolved) intervals are split as usual. Once the convergence is 2nd order, we are nearing the required resolution and can trust higher order methods. At this point, the trapezoid stencil is replaced with a Romberg tableau. When an Interval is processed, instead of splitting it's Romberg integrator is advanced one level to a maximum of 9.  This integrator is somewhat slower than the Gauss-Kronrod rules to get to machine precision, but offers robust error estimation for larger errors. In practice can perform faster and more robustly than other approaches if the required tolerance is only, say, `1.0e-3`.

The call signature is:

```c
double cadre_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, Mesh9 *mout, int verbose, int (*errf)(void *),
                  double *pfa, double *pfb)
```

The arguments are:

* `f`: A pointer to the function to integrate, called as `y = f(x, args)`.
* `xa`, `xb`: Left and right integration bounds
* `Nmax`: Maximum number of evaluations. 
* `atol`: Absolute error tolerance
* `rtol`: Relative error tolerance
* `args`: pointer to the arguments passed to `f`
* `errf`: pointer to the error function. Is called after every evaluation
    of `f`.  If return value is non-zero, integration is immediately halted
* `Neval`: pointer to integer. if not NULL, set to the number of evaluations performed
* `eps`: pointer to double. if not NULL, set to the error estimate
* `mout`: a pointer to a mesh or NULL. If not NULL, the entire integration mesh is returned. This is an allocated structure and must be freed with the appropriate method from `interval.h`. Either a `Mesh`, a `Mesh3`, a `Mesh5`, or a `Mesh9`.
* `verbose`: if non-zero, print internal information during run
* `errf`: pointer to the error function. Is called after every evaluation of `f`.  If return value is non-zero, integration is immediately halted.
* `pfa`, `pfb`: pointers to double. If not NULL, used for `f(xa)` and `f(xb)`.
