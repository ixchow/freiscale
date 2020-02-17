/******************************************************************************
*  DCT Consistency Check
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "otfft/otfft.h"
using namespace std;
using OTFFT::complex_t;
using OTFFT::simd_malloc;
using OTFFT::simd_free;

void print_err1(const char* name, int N, const double* x, const double* y)
{
    if (N < (1<<15)) {
        double* z  = (double*) simd_malloc(N*sizeof(double));
        for (int k = 0; k < N; k++) {
            z[k] = 0;
            for (int p = 0; p < N; p++)
                z[k] += x[p]*cos(M_PI/N*k*(p+1.0/2));
        }
        double err = 0;
        for (int k = 0; k < N; k++) {
            const double d = y[k] - z[k];
            err += d*d;
        }
        cout << name << ":";
        if (err == 0)
            cout << "---" << endl;
        else if (err > 0)
            cout << setw(3) << lrint(log10(err)) << endl;
        else
            cout << "NG" << endl;
        simd_free(z);
    }
}

void print_err2(const char* name, int N, const double* x, const double* y)
{
    double err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - x[k];
        err += d*d;
    }
    cout << name << ":";
    if (err == 0)
        cout << "---" << endl;
    else if (err > 0)
        cout << setw(3) << lrint(log10(err)) << endl;
    else
        cout << "NG" << endl;
}

void print_err3(const char* name, int N, const double* x, const double* y)
{
    if (N < (1<<15)) {
        double* z  = (double*) simd_malloc(N*sizeof(double));
        for (int k = 0; k < N; k++) {
            z[k] = 0;
            for (int p = 0; p < N; p++)
                z[k] += x[p]*cos(M_PI/N*k*(p+1.0/2));
            z[k] /= N;
        }
        double err = 0;
        for (int k = 0; k < N; k++) {
            const double d = y[k] - z[k];
            err += d*d;
        }
        cout << name << ":";
        if (err == 0)
            cout << "---" << endl;
        else if (err > 0)
            cout << setw(3) << lrint(log10(err)) << endl;
        else
            cout << "NG" << endl;
        simd_free(z);
    }
}

int main() try
{
    srand((unsigned) time(NULL));
    for (int n = 1; n <= 22; n++) {
        const int N = 1 << n;
        double* x0 = (double*) simd_malloc(N*sizeof(double));
        double* x  = (double*) simd_malloc(N*sizeof(double));
        for (int p = 0; p < N; p++) x0[p] = rand() % 100 - 50;
        OTFFT::DCT dct(N);
        cout << "[2^(" << setw(2) << n << ")]" << endl;

        for (int p = 0; p < N; p++) x[p] = x0[p];
        dct.fwd0(x);
        print_err1("fwd0", N, x0, x);
        dct.invn(x);
        print_err2("invn", N, x0, x);

        for (int p = 0; p < N; p++) x[p] = x0[p];
        dct.fwdn(x);
        print_err3("fwdn", N, x0, x);
        dct.inv0(x);
        print_err2("inv0", N, x0, x);

        for (int p = 0; p < N; p++) x[p] = x0[p];
        dct.fwd(x);
        print_err3("fwd ", N, x0, x);
        dct.inv(x);
        print_err2("inv ", N, x0, x);

        simd_free(x);
        simd_free(x0);
    }
    return 0;
}
catch (...) { cerr << "\n""*** exception! ***" << endl; }
