/******************************************************************************
*  Simple DFT
*
*  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_simpledft_h
#define otfft_simpledft_h

#include "otfft/otfft_misc.h"

namespace SimpleDFT { /////////////////////////////////////////////////////////

using namespace OTFFT_MISC;

///////////////////////////////////////////////////////////////////////////////

void dft(const int n,
        complex_vector x, complex_vector y, const_complex_vector W)
{
    typedef long long integer;
    for (integer k = 0; k < n; k++) {
        complex_t z = 0;
        for (integer p = 0; p < n; p++) z += x[p]*W[(k*p)%n];
        y[k] = z;
    }
    for (int k = 0; k < n; k++) x[k] = y[k];
}

void idft(const int n,
        complex_vector x, complex_vector y, const_complex_vector W)
{
    typedef long long integer;
    for (integer k = 0; k < n; k++) {
        complex_t z = 0;
        for (integer p = 0; p < n; p++) z += x[p]*W[n-(k*p)%n];
        y[k] = z;
    }
    for (int k = 0; k < n; k++) x[k] = y[k];
}

///////////////////////////////////////////////////////////////////////////////
// DFT object
///////////////////////////////////////////////////////////////////////////////

struct DFT0
{
    int N, log_N;
    simd_array<complex_t> weight;
    complex_t* __restrict W;

    DFT0() NOEXCEPT : N(0), W(0) {}
    DFT0(const int n) { setup(n); }

    void setup(const int n)
    {
        const double theta0 = 2*M_PI/n;
        N = n;
        weight.setup(n+1); W = &weight;
        for (int p = 0; p <= n; p++) {
            W[p] = complex_t(cos(p*theta0), -sin(p*theta0));
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    inline void fwd(complex_vector x, complex_vector y) const NOEXCEPT
    {
        dft(N, x, y, W);
        for (int k = 0; k < N; k++) x[k] /= N;
    }

    inline void fwd0(complex_vector x, complex_vector y) const NOEXCEPT
    {
        dft(N, x, y, W);
    }

    inline void fwdu(complex_vector x, complex_vector y) const NOEXCEPT
    {
        const double sN = sqrt(double(N));
        dft(N, x, y, W);
        for (int k = 0; k < N; k++) x[k] /= sN;
    }

    inline void fwdn(complex_vector x, complex_vector y) const NOEXCEPT
    {
        fwd(x, y);
    }

    ///////////////////////////////////////////////////////////////////////////

    inline void inv(complex_vector x, complex_vector y) const NOEXCEPT
    {
        idft(N, x, y, W);
    }

    inline void inv0(complex_vector x, complex_vector y) const NOEXCEPT
    {
        inv(x, y);
    }

    inline void invu(complex_vector x, complex_vector y) const NOEXCEPT
    {
        const double sN = sqrt(double(N));
        idft(N, x, y, W);
        for (int p = 0; p < N; p++) x[p] /= sN;
    }

    inline void invn(complex_vector x, complex_vector y) const NOEXCEPT
    {
        idft(N, x, y, W);
        for (int p = 0; p < N; p++) x[p] /= N;
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_simpledft_h
