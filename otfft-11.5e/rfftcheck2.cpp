/******************************************************************************
*  RFFT Consistency Check 2
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "simple_dft.h"
#include "otfft/otfft.h"
using namespace std;
using OTFFT::complex_t;
using OTFFT::simd_malloc;
using OTFFT::simd_free;

void print_err1(const char* name, int N, const complex_t* x, const complex_t* y)
{
    double err = 0;
    for (int k = 0; k < N; k++) {
        const complex_t d = y[k] - x[k];
        err += Re(d*conj(d));
    }
    cout << name << ":";
    if (err == 0)
        cout << "---" << endl;
    else if (err > 0)
        cout << setw(3) << lrint(log10(err)) << endl;
    else
        cout << "NG" << endl;
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

int power(const int a, const int n)
{
    int p = 1;
    for (int k = 0; k < n; k++) p *= a;
    return p;
}

int main() try
{
    srand((unsigned) time(NULL));
    for (int n = 1; n <= 4; n++) {
        const int N = power(10, n);
        //const int N = 2*power(3, n);
        double*    x0 = (double*) simd_malloc(N*sizeof(double));
        double*    x  = (double*) simd_malloc(N*sizeof(double));
        complex_t* y  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* z  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int p = 0; p < N; p++) x0[p] = rand() % 100 - 50;
        SimpleDFT::DFT0 simple_dft(N);
        OTFFT::RFFT otfft(N);
        cout << "[" << setw(6) << N << "]" << endl;

        for (int p = 0; p < N; p++) y[p] = x0[p];
        simple_dft.fwd0(y, z);
        otfft.fwd0(x0, z);
        print_err1("fwd0", N, y, z);
        otfft.invn(z, x);
        print_err2("invn", N, x0, x);

        for (int p = 0; p < N; p++) y[p] = x0[p];
        simple_dft.fwd(y, z);
        otfft.fwd(x0, z);
        print_err1("fwd ", N, y, z);
        otfft.inv(z, x);
        print_err2("inv ", N, x0, x);

        for (int p = 0; p < N; p++) y[p] = x0[p];
        simple_dft.fwdu(y, z);
        otfft.fwdu(x0, z);
        print_err1("fwdu", N, y, z);
        otfft.invu(z, x);
        print_err2("invu", N, x0, x);

        for (int p = 0; p < N; p++) y[p] = x0[p];
        simple_dft.fwdn(y, z);
        otfft.fwdn(x0, z);
        print_err1("fwdn", N, y, z);
        otfft.inv0(z, x);
        print_err2("inv0", N, x0, x);

        simd_free(z);
        simd_free(y);
        simd_free(x);
        simd_free(x0);
    }
    return 0;
}
catch (...) { cerr << "\n""*** exception! ***" << endl; }
