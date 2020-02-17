/******************************************************************************
*  FFT Consistency Check
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

void print_err1(const char* name, int N, const complex_t* x, const complex_t* y)
{
    if (N < (1<<15)) {
        complex_t* z  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int k = 0; k < N; k++) {
            z[k] = 0;
            for (int p = 0; p < N; p++) {
                const double theta = (2*M_PI/N)*k*p;
                z[k] += x[p]*complex_t(cos(theta), -sin(theta));
            }
        }
        double err = 0;
        for (int k = 0; k < N; k++) {
            const complex_t d = y[k] - z[k];
            err += Re(d*conj(d));
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

void print_err2(const char* name, int N, const complex_t* x, const complex_t* y)
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

void print_err3(const char* name, int N, const complex_t* x, const complex_t* y)
{
    if (N < (1<<15)) {
        complex_t* z  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int k = 0; k < N; k++) {
            z[k] = 0;
            for (int p = 0; p < N; p++) {
                const double theta = (2*M_PI/N)*k*p;
                z[k] += x[p]*complex_t(cos(theta), -sin(theta));
            }
            z[k] /= N;
        }
        double err = 0;
        for (int k = 0; k < N; k++) {
            const complex_t d = y[k] - z[k];
            err += Re(d*conj(d));
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

void print_err4(const char* name, int N, const complex_t* x, const complex_t* y)
{
    if (N < (1<<15)) {
        complex_t* z  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int k = 0; k < N; k++) {
            z[k] = 0;
            for (int p = 0; p < N; p++) {
                const double theta = (2*M_PI/N)*k*p;
                z[k] += x[p]*complex_t(cos(theta), -sin(theta));
            }
            z[k] /= sqrt(N);
        }
        double err = 0;
        for (int k = 0; k < N; k++) {
            const complex_t d = y[k] - z[k];
            err += Re(d*conj(d));
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
        const int N = (1 << n) + 123;
        complex_t* x0 = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* x  = (complex_t*) simd_malloc(N*sizeof(complex_t));
        for (int p = 0; p < N; p++) {
            x0[p].Re = rand() % 100 - 50;
            x0[p].Im = rand() % 100 - 50;
        }
        OTFFT::Bluestein bst(N);
        cout << "[" << setw(7) << N << "]" << endl;

        for (int p = 0; p < N; p++) x[p] = x0[p];
        bst.fwd0(x);
        print_err1("fwd0", N, x0, x);
        bst.invn(x);
        print_err2("invn", N, x0, x);
        for (int p = 0; p < N; p++) x[p] = x0[p];
        bst.fwdn(x);
        print_err3("fwdn", N, x0, x);
        bst.inv0(x);
        print_err2("inv0", N, x0, x);
        for (int p = 0; p < N; p++) x[p] = x0[p];
        bst.fwd(x);
        print_err3("fwd ", N, x0, x);
        bst.inv(x);
        print_err2("inv ", N, x0, x);
        for (int p = 0; p < N; p++) x[p] = x0[p];
        bst.fwdu(x);
        print_err4("fwdu", N, x0, x);
        bst.invu(x);
        print_err2("invu", N, x0, x);

        simd_free(x);
        simd_free(x0);
    }
    return 0;
}
catch (...) { cerr << "\n""*** exception! ***" << endl; }
