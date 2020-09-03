#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integrate.h"
#include "interval.h"


double f1(double x, void *args)
{
    double c = ((double *)args)[0];

    return 1.0 / (1.0 + x*x/(c*c));
}

double F1(double a, double b, void *args)
{
    double c = ((double *)args)[0];

    return c * (atan(b/c) - atan(a/c));
}

double f2(double x, void *args)
{
    double c = ((double *)args)[0];

    return cos(c*x);
}

double F2(double a, double b, void *args)
{
    double c = ((double *)args)[0];

    return (sin(c*b) - sin(c*a)) / c;
}

void dumpConvergence(char name[], int N, int *Neval, double *I, double *err)
{
    char fname[256];
    sprintf(fname, "convergence_%s.txt", name);

    FILE *f = fopen(fname, "w");
    fprintf(f, "%s\n", name);

    int i;
    for(i=0; i<N; i++)
        fprintf(f, "%d %.16lg %.16lg\n", Neval[i], I[i], err[i]);

    fclose(f);
}

void dumpConvergenceAdapt(char name[], int N, int *Neval, double *I,
                          double *err, double *err_adapt)
{
    char fname[256];
    sprintf(fname, "convergence_%s.txt", name);

    FILE *f = fopen(fname, "w");
    fprintf(f, "%s\n", name);

    int i;
    for(i=0; i<N; i++)
        fprintf(f, "%d %.16lg %.16lg %.16lg\n", Neval[i], I[i],
                err[i], err_adapt[i]);

    fclose(f);
}

void convergence_trap_fixed(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        Neval[i] = n+1;
        I[i] = trap(f, a, b, n, args);
        err[i] = fabs((I[i] - exact)/exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergence(name, N, Neval, I, err);
}

void convergence_simp_fixed(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        if(n%2 == 1)
            n--;

        Neval[i] = n+1;
        I[i] = simp(f, a, b, n, args);
        err[i] = fabs((I[i] - exact)/exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergence(name, N, Neval, I, err);
}

void convergence_romb_adapt(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);
    double err_approx[N];

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        I[i] = romb(f, a, b, n, 0, 0, args, &(Neval[i]), &(err_approx[i]), 0);
        err[i] = fabs((I[i] - exact)/exact);
        err_approx[i] = fabs(err_approx[i] / exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergenceAdapt(name, N, Neval, I, err, err_approx);
}

void convergence_trap_adapt(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);
    double err_approx[N];

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        I[i] = trap_adapt(f, a, b, n, 0, 0, args,
                          &(Neval[i]),  &(err_approx[i]), NULL, 0);
        err[i] = fabs((I[i] - exact)/exact);
        err_approx[i] = fabs(err_approx[i] / exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergenceAdapt(name, N, Neval, I, err, err_approx);
}

void convergence_simp_adapt(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);

    double err_approx[N];

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        I[i] = simp_adapt(f, a, b, n, 0, 0, args,
                          &(Neval[i]), &(err_approx[i]), NULL, 0);
        err[i] = fabs((I[i] - exact)/exact);
        err_approx[i] = fabs(err_approx[i] / exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergenceAdapt(name, N, Neval, I, err, err_approx);
}

void convergence_gl7_adapt(double (*f)(double, void *), void *args,
                            double a, double b,
                            double (*F)(double, double, void *),
                            int N1, int N2,
                            int N, int *Neval, double *I, double *err,
                            char name[])
{
    double exact = F(a, b, args);

    double err_approx[N];

    int i;
    for(i=0; i<N; i++)
    {
        int n = (int) (N1 * pow(((double)N2)/N1, ((double) i)/(N-1)));
        if(i==0)
            n = N1;
        else if(i == N-1)
            n = N2;

        I[i] = gl7_adapt(f, a, b, n, 0, 0, args, &(Neval[i]), &(err_approx[i]),
                         NULL, 0);
        err[i] = fabs((I[i] - exact)/exact);
        err_approx[i] = fabs(err_approx[i] / exact);
        printf("%05d: %.12lg %.6e\n", Neval[i], I[i], err[i]);
    }

    dumpConvergenceAdapt(name, N, Neval, I, err, err_approx);
}

int main(int argc, char *argv[])
{
    double (*f)(double, void *);
    double (*F)(double, double, void *);
    double a = 0.0;
    double b = 1.0;
    double c = 0.0001;

    double args[1] = {c};

    f = f1;
    F = F1;

    double exact = F(a, b, args);

    int i;
    for(i=1; i<=7; i++)
    {
        double I = gl(f, a, b, i, args);
        printf("%d: %.12lf %.6le\n", i, I, (I-exact)/exact);
    }
    printf("   %.12lf\n", exact);

    int Neval;
    double err;
    struct Mesh3 m;
    double I = trap_adapt(f, a, b, 200, 1.0e-10, 1.0e-10, args,
                          &Neval, &err, &m, 1);
    double true_error = I - exact;
    printf("   %.12lf %.6le %.6le %d\n", I, err, true_error, Neval);
    printf("%lu %lu\n", m.N, m.totalSize);
    /*
    for(i=0; i<m.N; i++)
        printf("%d: (%.6lf, %.6lf) %.6lf %.3le\n", i, m.heap[i].a, m.heap[i].b,
               m.heap[i].I, m.heap[i].err);
    */
    mesh3Free(&m);
    
    
    struct Mesh m2;
    I = gl7_adapt(f, a, b, 200, 1.0e-10, 1.0e-10, args,
                  &Neval, &err, &m2, 1);
    true_error = I - exact;
    printf("   %.12lf %.6le %.6le %d\n", I, err, true_error, Neval);
    printf("%lu %lu\n", m2.N, m2.totalSize);
    /*
    for(i=0; i<m2.N; i++)
        printf("%d: (%.6lf, %.6lf) %.6lf %.3le\n", i,
               m2.heap[i].a, m2.heap[i].b,
               m2.heap[i].I, m2.heap[i].err);
               */
    meshFree(&m2);

    int N = 500;
    int tf_Neval[N];
    double tf_I[N];
    double tf_err[N];
    int ta_Neval[N];
    double ta_I[N];
    double ta_err[N];
    int gl7a_Neval[N];
    double gl7a_I[N];
    double gl7a_err[N];

    convergence_trap_fixed(f, args, a, b, F, 2, 4096,
                           N, tf_Neval, tf_I, tf_err, "trap_fixed");
    convergence_trap_adapt(f, args, a, b, F, 2, 4096,
                           N, ta_Neval, ta_I, ta_err, "trap_adapt");
    convergence_simp_fixed(f, args, a, b, F, 2, 4096,
                           N, tf_Neval, tf_I, tf_err, "simp_fixed");
    convergence_simp_adapt(f, args, a, b, F, 2, 4096,
                           N, ta_Neval, ta_I, ta_err, "simp_adapt");
    convergence_romb_adapt(f, args, a, b, F, 2, 4097,
                           N, ta_Neval, ta_I, ta_err, "romb_adapt");
    convergence_gl7_adapt(f, args, a, b, F, 2, 4096,
                          N, gl7a_Neval, gl7a_I, gl7a_err, "gl7_adapt");

    return 0;
}
