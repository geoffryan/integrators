#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "integrate.h"
#include "interval.h"

#define KMAX 20

double trap(double (*f)(double, void *), double xa, double xb, int N, void *args)
{
    double dx = (xb - xa)/N;
    double I = 0.5*(f(xa, args) + f(xb, args));
    int i;
    for(i=1; i<N; i++)
        I += f(xa + i*dx, args);
    return I*dx;
}

double simp(double (*f)(double, void *), double xa, double xb, int N, void *args)
{
    if(N%2 == 1)
        N -= 1;

    double dx = (xb - xa)/N;
    double I1, I2, I3;
    int i;
    I1 = f(xa, args) + f(xb, args);
    I2 = 0.0;
    for(i=1; i<N; i+=2)
        I2 += f(xa + i*dx, args);
    I3 = 0.0;
    for(i=2; i<N; i+=2)
        I3 += f(xa + i*dx, args);
    return (I1 + 4*I2 + 2*I3) * dx / 3.0;
}

double romb(double (*f)(double, void *), double xa, double xb, int N,
            double atol, double rtol, void *args, int *Neval, double *eps,
            int verbose)
{
    double R[KMAX];

    int m, k, k0, Nk;
    long fpm;
    double hk, Rp, err;

    double maxFracChange = 0.1;

    hk = xb - xa;
    Nk = 1;
    R[KMAX-1] = 0.5*(xb-xa)*(f(xa, args) + f(xb, args));
    R[0] = R[KMAX-1];

    if(Neval != NULL)
        *Neval = 2;

    for(k=1; k<KMAX; k++)
    {
        k0 = KMAX-k-1;
        hk *= 0.5;
        Nk *= 2;

        Rp = 0.0;
        for(m=1; m<Nk; m+=2)
            Rp += f(xa + m*hk, args);
        R[k0] = 0.5*R[k0+1] + hk*Rp;
        if(Neval != NULL)
            *Neval += Nk/2;

        fpm = 1;
        for(m=1; m<=k; m++)
        {
            fpm *= 4;
            R[k0+m] = (fpm*R[k0+m-1] - R[k0+m]) / (fpm - 1);
        }
        err = (R[KMAX-1] - R[0]) / (fpm - 1);
        double lastVal = R[0];
        R[0] = R[KMAX-1];

        double fracChange = fabs((R[0] - lastVal) / lastVal);

        if(eps != NULL)
            *eps = err;
        if(verbose)
            printf("level %d:  Neval=%d  I=%.6lg  fracDelta=%.3lg"
                   " err=%.6lg  tol=%.6lg\n",
                    k, Nk+1, R[0], fracChange, err, atol+rtol*fabs(R[0]));

        //printf("      k%d: I=%.6le err=%.6le frac=%.3le\n", k, R[0], err, 
        //        fabs(err) / (atol + rtol*fabs(R[0])));

        if((fabs(err) < atol + rtol*fabs(R[0]))
                && fracChange < maxFracChange)
            break;

        if(N > 1 && Nk >= N)
            break;
    }

    return R[0];
}


double gl(double (*f)(double, void *), double xa, double xb, int N,
          void *args)
{
    double zi[N];
    double wi[N];
    glPoints(&N, zi, wi);
    return quad(f, xa, xb, N, args, zi, wi, -1, 1);
}

double quad(double (*f)(double, void *), double xa, double xb, int N,
            void *args, double *zi, double *wi, double za, double zb)
{
    double m = (xb - xa) / (zb - za);
    double c = xa - m*za;

    int i;

    double I = 0;
    for(i=0; i<N; i++)
        I += wi[i] * f(m*zi[i]+c, args);

    I *= m;

    return I;

}

double trap_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, Mesh3 *mout, int verbose)
{
    double I = m3_adapt(f, xa, xb, Nmax, trapProcessInterval,
                        trapSplitInterval, atol, rtol, args,
                        Neval, eps, mout, verbose);
    return I;
}

double gl7_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh *mout, int verbose)
{
    double I = gen_adapt(f, xa, xb, Nmax, gl7ProcessInterval, gl7SplitInterval,
                         atol, rtol, args, Neval, eps, mout, verbose);
    return I;
}

double gen_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval *),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                      Interval *, Interval *, Interval *),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh *mout, int verbose)
{
    Mesh m;
    meshInit(&m);

    Interval i = {.a=xa, .b=xb, .I=0, .err=0};
    int n = processInterval(f, args, &i);

    meshInsert(&m, i);

    double I = i.I;
    double err = i.err;
    int num_intervals = 1;

    while(n < Nmax)
    {
        i = meshExtract(&m);

        Interval i1;
        Interval i2;
        n += splitInterval(f, args, &i, &i1, &i2);
        meshInsert(&m, i1);
        meshInsert(&m, i2);
        num_intervals++;


        err += i1.err + i2.err - i.err;
        I += i1.I + i2.I - i.I;
        
        if(verbose)
            printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                   "  meshOk=%d\n",
                   num_intervals, I, err, atol + rtol*fabs(I), meshCheck(&m));

        if(err < atol + rtol*fabs(I))
            break;
    }

    I = meshTotalIntegral(&m);

    if(Neval != NULL)
        *Neval = n;

    if(eps != NULL)
    {
        err = meshTotalError(&m);
        *eps = err;
    }

    if(mout == NULL)
        meshFree(&m);
    else
        *mout = m;

    return I;
}

double m3_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval3 *),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                      Interval3 *, Interval3 *, Interval3 *),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh3 *mout, int verbose)
{
    Mesh3 m;
    mesh3Init(&m);

    Interval3 i = {.a=xa, .b=xb, .I=0, .err=0, .fa=0, .fm=0, .fb=0};
    i.fa = f(xa, args);
    i.fb = f(xb, args);
    int n = 2;
    n += processInterval(f, args, &i);

    mesh3Insert(&m, i);

    double I = i.I;
    double err = i.err;
    int num_intervals = 1;

    while(n < Nmax)
    {
        i = mesh3Extract(&m);

        Interval3 i1;
        Interval3 i2;
        n += splitInterval(f, args, &i, &i1, &i2);
        mesh3Insert(&m, i1);
        mesh3Insert(&m, i2);
        num_intervals++;


        err += i1.err + i2.err - i.err;
        I += i1.I + i2.I - i.I;
        
        if(verbose)
            printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                   "  meshOk=%d\n",
                   num_intervals, I, err, atol + rtol*fabs(I), mesh3Check(&m));

        if(err < atol + rtol*fabs(I))
            break;
    }

    I = mesh3TotalIntegral(&m);

    if(Neval != NULL)
        *Neval = n;

    if(eps != NULL)
    {
        err = mesh3TotalError(&m);
        *eps = err;
    }

    if(mout == NULL)
        mesh3Free(&m);
    else
        *mout = m;

    return I;
}

double m5_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval5 *),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                      Interval5 *, Interval5 *, Interval5 *),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh5 *mout, int verbose)
{
    Mesh5 m;
    mesh5Init(&m);

    Interval5 i = {.a=xa, .b=xb, .I=0, .err=0,
                   .fa=0, .fl=0, .fm=0, .fr=0, .fb=0};
    int n = processInterval(f, args, &i);

    mesh5Insert(&m, i);

    double I = i.I;
    double err = i.err;
    int num_intervals = 1;

    while(n < Nmax)
    {
        i = mesh5Extract(&m);

        Interval5 i1;
        Interval5 i2;
        n += splitInterval(f, args, &i, &i1, &i2);
        mesh5Insert(&m, i1);
        mesh5Insert(&m, i2);
        num_intervals++;


        err += i1.err + i2.err - i.err;
        I += i1.I + i2.I - i.I;
        
        if(verbose)
            printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                   "  meshOk=%d\n",
                   num_intervals, I, err, atol + rtol*fabs(I), mesh5Check(&m));

        if(err < atol + rtol*fabs(I))
            break;
    }

    I = mesh5TotalIntegral(&m);

    if(Neval != NULL)
        *Neval = n;

    if(eps != NULL)
    {
        err = mesh5TotalError(&m);
        *eps = err;
    }

    if(mout == NULL)
        mesh5Free(&m);
    else
        *mout = m;

    return I;
}

int trapProcessInterval(double (*f)(double, void *), void *args, Interval3 *i)
{
    double fa = i->fa;
    double fb = i->fb;
    double fm = f(0.5*(i->a+i->b), args);
    i->fm = fm;

    double h = 0.5*(i->b - i->a);

    double I0 = 0.5*(fa + fb) * 2*h;
    double I1 = 0.5*(fa + 2*fm + fb) * h;

    double err = (I1 - I0) / 3.0;
    i->err = fabs(err);
    i->I = I1 + err;

    return 1;
}

int trapSplitInterval(double (*f)(double, void *), void *args,
                      Interval3 *i0, Interval3 *i1, Interval3 *i2)
{
    double x = 0.5*(i0->a + i0->b);
    i1->a = i0->a;
    i1->b = x;
    i2->a = x;
    i2->b = i0->b;

    i1->fa = i0->fa;
    i1->fb = i0->fm;
    i2->fa = i0->fm;
    i2->fb = i0->fb;

    int n = 0;
    n += trapProcessInterval(f, args, i1);
    n += trapProcessInterval(f, args, i2);

    return n;
}

int gl7ProcessInterval(double (*f)(double, void *), void *args, Interval *i)
{
    int N = 7;
    double z[7];
    double w[7];
    glPoints(&N, z, w);

    double I0 = quad(f, i->a, i->b, N, args, z, w, -1, 1);
    
    double x = 0.5*(i->a+i->b);
    double I1 = quad(f, i->a, x, N, args, z, w, -1, 1);
    double I2 = quad(f, x, i->b, N, args, z, w, -1, 1);

    i->err = fabs(I1 + I2 - I0);
    i->I = I1 + I2;

    return 3*N;
}

int gl7SplitInterval(double (*f)(double, void *), void *args,
                     Interval *i0, Interval *i1, Interval *i2)
{
    double x = 0.5*(i0->a + i0->b);
    i1->a = i0->a;
    i1->b = x;
    i2->a = x;
    i2->b = i0->b;

    int n = 0;
    n += gl7ProcessInterval(f, args, i1);
    n += gl7ProcessInterval(f, args, i2);

    return n;
}

/*
double adapt_general()
{
    intervalHeap[];

    err = total_err(intervalHeap);

    while(err > tol)
    {
        worstInterval = popMaxErr(intervalHeap);
        errWorst = worstInterval.error
        intervalA, intervalB = subdivide(worstInterval)
        process(intervalA)
        process(intervalB)
        intervalHeap.add(intervalA, intervalB)
        errA, errB = intervalA.error, intervalB.error
        err += (errA+errB-errWorst)
    }
}
*/

void glPoints(int *N, double *x, double *w)
{
    if(*N < 1)
        return;
    else if(*N == 1)
    {
        x[0] = 0.0;
        w[0] = 2.0;
    }
    else if(*N == 2)
    {
        x[0] = -0.5773502691896257;
        x[1] = 0.5773502691896257;
        w[0] = 1.0;
        w[1] = 1.0;
    }
    else if(*N == 3)
    {
        x[0] = -0.7745966692414834;
        x[1] = 0.0;
        x[2] = 0.7745966692414834;
        w[0] = 0.5555555555555556;
        w[1] = 0.8888888888888888;
        w[2] = 0.5555555555555556;
    }
    else if(*N == 4)
    {
        x[0] = -0.8611363115940526;
        x[1] = -0.3399810435848563;
        x[2] = 0.3399810435848563;
        x[3] = 0.8611363115940526;
        w[0] = 0.3478548451374538;
        w[1] = 0.6521451548625461;
        w[2] = 0.6521451548625461;
        w[3] = 0.3478548451374538;
    }
    else if(*N == 5)
    {
        x[0] = -0.9061798459386640;
        x[1] = -0.5384693101056831;
        x[2] = 0.0;
        x[3] = 0.5384693101056831;
        x[4] = 0.9061798459386640;
        w[0] = 0.2369268850561891;
        w[1] = 0.4786286704993665;
        w[2] = 0.5688888888888889;
        w[3] = 0.4786286704993665;
        w[4] = 0.2369268850561891;
    }
    else if(*N == 6)
    {
        x[0] = -0.9324695142031521;
        x[1] = -0.6612093864662645;
        x[2] = -0.2386191860831969;
        x[3] = 0.2386191860831969;
        x[4] = 0.6612093864662645;
        x[5] = 0.9324695142031521;
        w[0] = 0.1713244923791704;
        w[1] = 0.3607615730481386;
        w[2] = 0.4679139345726910;
        w[3] = 0.4679139345726910;
        w[4] = 0.3607615730481386;
        w[5] = 0.1713244923791704;
    }
    else if(*N == 7)
    {
        x[0] = -0.9491079123427585;
        x[1] = -0.7415311855993945;
        x[2] = -0.4058451513773972;
        x[3] = 0.0;
        x[4] = 0.4058451513773972;
        x[5] = 0.7415311855993945;
        x[6] = 0.9491079123427585;
        w[0] = 0.1294849661688697;
        w[1] = 0.2797053914892766;
        w[2] = 0.3818300505051189;
        w[3] = 0.4179591836734694;
        w[4] = 0.3818300505051189;
        w[5] = 0.2797053914892766;
        w[6] = 0.1294849661688697;
    }
    else
    {
        x[0] = -0.9491079123427585;
        x[1] = -0.7415311855993945;
        x[2] = -0.4058451513773972;
        x[3] = 0.0;
        x[4] = 0.4058451513773972;
        x[5] = 0.7415311855993945;
        x[6] = 0.9491079123427585;
        w[0] = 0.1294849661688697;
        w[1] = 0.2797053914892766;
        w[2] = 0.3818300505051189;
        w[3] = 0.4179591836734694;
        w[4] = 0.3818300505051189;
        w[5] = 0.2797053914892766;
        w[6] = 0.1294849661688697;
        *N = 7;
    }
}
