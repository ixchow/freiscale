/******************************************************************************
*  FFT Consistency Check 2
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

int power(const int a, const int n)
{
    int p = 1;
    for (int k = 0; k < n; k++) p *= a;
    return p;
}

void print_err(const char* name, int N, const complex_t* x, const complex_t* y)
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

int main() try
{
    srand((unsigned) time(NULL));
    for (int n = 1; n <= 4; n++) {
        const int N = power(10, n);
        //const int N = 3*power(7, n);
        //const int N = 2*3*power(7, n);
        complex_t* x0 = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* x  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* y  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* z  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int p = 0; p < N; p++) {
            x0[p].Re = rand() % 100 - 50;
            x0[p].Im = rand() % 100 - 50;
        }
        SimpleDFT::DFT0 simple_dft(N);
        OTFFT::FFT otfft(N);
        cout << "[" << setw(7) << N << "]" << endl;

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.fwd(x, z);
        otfft.fwd(y);
        print_err("fwd ", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.fwd0(x, z);
        otfft.fwd0(y);
        print_err("fwd0", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.fwdu(x, z);
        otfft.fwdu(y);
        print_err("fwdu", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.fwdn(x, z);
        otfft.fwdn(y);
        print_err("fwdn", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.inv(x, z);
        otfft.inv(y);
        print_err("inv ", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.inv0(x, z);
        otfft.inv0(y);
        print_err("inv0", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.invu(x, z);
        otfft.invu(y);
        print_err("invu", N, x, y);

        for (int p = 0; p < N; p++) x[p] = y[p] = x0[p];
        simple_dft.invn(x, z);
        otfft.invn(y);
        print_err("invn", N, x, y);

        simd_free(y);
        simd_free(x);
        simd_free(x0);
    }
    return 0;
}
catch (...) { cerr << "\n""*** exception! ***" << endl; }
