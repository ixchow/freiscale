/******************************************************************************
*  Simple FFT 2(Cooley Tukey Radix-4)
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef simple_fft2_h
#define simple_fft2_h

#include <utility>
#include "otfft/otfft_misc.h"

namespace SimpleFFT { /////////////////////////////////////////////////////////

using namespace OTFFT_MISC;

#ifdef DO_SINGLE_THREAD
static const int OMP_THRESHOLD = 1<<30;
#else
static const int OMP_THRESHOLD = 1<<15;
#endif

void fwdbut(int N, complex_vector x, const_complex_vector W) NOEXCEPT
{
    int n = N;
    for (int s = 1; n > 2; n /= 4, s *= 4) {
        const int n0 = 0;
        const int n1 = n / 4;
        const int n2 = n / 2;
        const int n3 = n1 + n2;
        if (N < OMP_THRESHOLD) {
            for (int q = 0; q < N; q += n) {
                complex_vector xq = x + q;
                if (n1 >= 2) for (int p = 0; p < n1; p += 2) {
                    const int sp = s*p;
                    const ymm w1p = cmplx2(W[sp], W[sp+s]);
                    const ymm w2p = mulpz2(w1p, w1p);
                    const ymm w3p = mulpz2(w1p, w2p);
                    const ymm a = getpz2(xq+p+n0);
                    const ymm b = getpz2(xq+p+n1);
                    const ymm c = getpz2(xq+p+n2);
                    const ymm d = getpz2(xq+p+n3);
                    const ymm  apc =       addpz2(a, c);
                    const ymm  amc =       subpz2(a, c);
                    const ymm  bpd =       addpz2(b, d);
                    const ymm jbmd = jxpz2(subpz2(b, d));
                    setpz2(xq+p+n0,             addpz2(apc,  bpd));
                    setpz2(xq+p+n1, mulpz2(w2p, subpz2(apc,  bpd)));
                    setpz2(xq+p+n2, mulpz2(w1p, subpz2(amc, jbmd)));
                    setpz2(xq+p+n3, mulpz2(w3p, addpz2(amc, jbmd)));
                }
                else {
                    const xmm a = getpz(xq[0]);
                    const xmm b = getpz(xq[1]);
                    const xmm c = getpz(xq[2]);
                    const xmm d = getpz(xq[3]);
                    const xmm  apc =      addpz(a, c);
                    const xmm  amc =      subpz(a, c);
                    const xmm  bpd =      addpz(b, d);
                    const xmm jbmd = jxpz(subpz(b, d));
                    setpz(xq[0], addpz(apc,  bpd));
                    setpz(xq[1], subpz(apc,  bpd));
                    setpz(xq[2], subpz(amc, jbmd));
                    setpz(xq[3], addpz(amc, jbmd));
                }
            }
        }
        else {
            if (n1 >= 2) {
                #pragma omp parallel for schedule(static)
                for (int q = 0; q < N; q += n) {
                    complex_vector xq = x + q;
                    for (int p = 0; p < n1; p += 2) {
                        const int sp = s*p;
                        const ymm w1p = cmplx2(W[sp], W[sp+s]);
                        const ymm w2p = mulpz2(w1p, w1p);
                        const ymm w3p = mulpz2(w1p, w2p);
                        const ymm a = getpz2(xq+p+n0);
                        const ymm b = getpz2(xq+p+n1);
                        const ymm c = getpz2(xq+p+n2);
                        const ymm d = getpz2(xq+p+n3);
                        const ymm  apc =       addpz2(a, c);
                        const ymm  amc =       subpz2(a, c);
                        const ymm  bpd =       addpz2(b, d);
                        const ymm jbmd = jxpz2(subpz2(b, d));
                        setpz2(xq+p+n0,             addpz2(apc,  bpd));
                        setpz2(xq+p+n1, mulpz2(w2p, subpz2(apc,  bpd)));
                        setpz2(xq+p+n2, mulpz2(w1p, subpz2(amc, jbmd)));
                        setpz2(xq+p+n3, mulpz2(w3p, addpz2(amc, jbmd)));
                    }
                }
            }
            else {
                #pragma omp parallel for schedule(static)
                for (int q = 0; q < N; q += n) {
                    complex_vector xq = x + q;
                    const xmm a = getpz(xq[0]);
                    const xmm b = getpz(xq[1]);
                    const xmm c = getpz(xq[2]);
                    const xmm d = getpz(xq[3]);
                    const xmm  apc =      addpz(a, c);
                    const xmm  amc =      subpz(a, c);
                    const xmm  bpd =      addpz(b, d);
                    const xmm jbmd = jxpz(subpz(b, d));
                    setpz(xq[0], addpz(apc,  bpd));
                    setpz(xq[1], subpz(apc,  bpd));
                    setpz(xq[2], subpz(amc, jbmd));
                    setpz(xq[3], addpz(amc, jbmd));
                }
            }
        }
    }
    if (n == 2) {
        if (N < OMP_THRESHOLD) {
            for (int q = 0; q < N; q += 2) {
                complex_t* const xq = x + q;
                const xmm a = getpz(xq[0]);
                const xmm b = getpz(xq[1]);
                setpz(xq[0], addpz(a, b));
                setpz(xq[1], subpz(a, b));
            }
        }
        else {
            #pragma omp parallel for schedule(static)
            for (int q = 0; q < N; q += 2) {
                complex_t* const xq = x + q;
                const xmm a = getpz(xq[0]);
                const xmm b = getpz(xq[1]);
                setpz(xq[0], addpz(a, b));
                setpz(xq[1], subpz(a, b));
            }
        }
    }
}

void invbut(int N, complex_vector x, const_complex_vector W) NOEXCEPT
{
    int n = N;
    for (int s = 1; n > 2; n /= 4, s *= 4) {
        const int n0 = 0;
        const int n1 = n / 4;
        const int n2 = n / 2;
        const int n3 = n1 + n2;
        if (N < OMP_THRESHOLD) {
            for (int q = 0; q < N; q += n) {
                complex_vector xq = x + q;
                if (n1 >= 2) for (int p = 0; p < n1; p += 2) {
                    const int sp = s*p;
                    const ymm w1p = cmplx2(W[N-sp], W[N-sp-s]);
                    const ymm w2p = mulpz2(w1p, w1p);
                    const ymm w3p = mulpz2(w1p, w2p);
                    const ymm a = getpz2(xq+p+n0);
                    const ymm b = getpz2(xq+p+n1);
                    const ymm c = getpz2(xq+p+n2);
                    const ymm d = getpz2(xq+p+n3);
                    const ymm  apc =       addpz2(a, c);
                    const ymm  amc =       subpz2(a, c);
                    const ymm  bpd =       addpz2(b, d);
                    const ymm jbmd = jxpz2(subpz2(b, d));
                    setpz2(xq+p+n0,             addpz2(apc,  bpd));
                    setpz2(xq+p+n1, mulpz2(w2p, subpz2(apc,  bpd)));
                    setpz2(xq+p+n2, mulpz2(w1p, addpz2(amc, jbmd)));
                    setpz2(xq+p+n3, mulpz2(w3p, subpz2(amc, jbmd)));
                }
                else {
                    const xmm a = getpz(xq[0]);
                    const xmm b = getpz(xq[1]);
                    const xmm c = getpz(xq[2]);
                    const xmm d = getpz(xq[3]);
                    const xmm  apc =      addpz(a, c);
                    const xmm  amc =      subpz(a, c);
                    const xmm  bpd =      addpz(b, d);
                    const xmm jbmd = jxpz(subpz(b, d));
                    setpz(xq[0], addpz(apc,  bpd));
                    setpz(xq[1], subpz(apc,  bpd));
                    setpz(xq[2], addpz(amc, jbmd));
                    setpz(xq[3], subpz(amc, jbmd));
                }
            }
        }
        else {
            if (n1 >= 2) {
                #pragma omp parallel for schedule(static)
                for (int q = 0; q < N; q += n) {
                    complex_vector xq = x + q;
                    for (int p = 0; p < n1; p += 2) {
                        const int sp = s*p;
                        const ymm w1p = cmplx2(W[N-sp], W[N-sp-s]);
                        const ymm w2p = mulpz2(w1p, w1p);
                        const ymm w3p = mulpz2(w1p, w2p);
                        const ymm a = getpz2(xq+p+n0);
                        const ymm b = getpz2(xq+p+n1);
                        const ymm c = getpz2(xq+p+n2);
                        const ymm d = getpz2(xq+p+n3);
                        const ymm  apc =       addpz2(a, c);
                        const ymm  amc =       subpz2(a, c);
                        const ymm  bpd =       addpz2(b, d);
                        const ymm jbmd = jxpz2(subpz2(b, d));
                        setpz2(xq+p+n0,             addpz2(apc,  bpd));
                        setpz2(xq+p+n1, mulpz2(w2p, subpz2(apc,  bpd)));
                        setpz2(xq+p+n2, mulpz2(w1p, addpz2(amc, jbmd)));
                        setpz2(xq+p+n3, mulpz2(w3p, subpz2(amc, jbmd)));
                    }
                }
            }
            else {
                #pragma omp parallel for schedule(static)
                for (int q = 0; q < N; q += n) {
                    complex_vector xq = x + q;
                    const xmm a = getpz(xq[0]);
                    const xmm b = getpz(xq[1]);
                    const xmm c = getpz(xq[2]);
                    const xmm d = getpz(xq[3]);
                    const xmm  apc =      addpz(a, c);
                    const xmm  amc =      subpz(a, c);
                    const xmm  bpd =      addpz(b, d);
                    const xmm jbmd = jxpz(subpz(b, d));
                    setpz(xq[0], addpz(apc,  bpd));
                    setpz(xq[1], subpz(apc,  bpd));
                    setpz(xq[2], addpz(amc, jbmd));
                    setpz(xq[3], subpz(amc, jbmd));
                }
            }
        }
    }
    if (n == 2) {
        if (N < OMP_THRESHOLD) {
            for (int q = 0; q < N; q += 2) {
                complex_t* const xq = x + q;
                const xmm a = getpz(xq[0]);
                const xmm b = getpz(xq[1]);
                setpz(xq[0], addpz(a, b));
                setpz(xq[1], subpz(a, b));
            }
        }
        else {
            #pragma omp parallel for schedule(static)
            for (int q = 0; q < N; q += 2) {
                complex_t* const xq = x + q;
                const xmm a = getpz(xq[0]);
                const xmm b = getpz(xq[1]);
                setpz(xq[0], addpz(a, b));
                setpz(xq[1], subpz(a, b));
            }
        }
    }
}

struct FFT
{
    const int N;
    simd_array<complex_t> weight;
    simd_array<int> table;
    complex_t* const W;
    int* bitrev;

    FFT(int n) : N(n), weight(n+1), table(n), W(&weight), bitrev(&table)
    {
        init_W(N, W);
        bitrev[0] = 0; bitrev[N-1] = N-1;
        for (int i = 0, j = 1; j < N-1; j++) {
            for (int k = N >> 1; k > (i ^= k); k >>= 1);
            bitrev[j] = i;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    void fwd(complex_vector x) const NOEXCEPT
    {
        const ymm rN = { 1.0/N, 1.0/N, 1.0/N, 1.0/N };
        fwdbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            if (N >= 2) {
                for (int p = 0; p < N; p += 2) {
                    setpz2(x+p, mulpd2(rN, getpz2(x+p)));
                }
            }
        }
        else
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N; p += 2) {
                setpz2(x+p, mulpd2(rN, getpz2(x+p)));
            }
        }
    }

    void fwd0(complex_vector x) const NOEXCEPT
    {
        fwdbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
        }
        else {
            #pragma omp parallel for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
        }
    }

    void fwdu(complex_vector x) const NOEXCEPT
    {
        const double ssrN = sqrt(1.0/N);
        const ymm srN = { ssrN, ssrN, ssrN, ssrN };
        fwdbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            if (N >= 2) {
                for (int p = 0; p < N; p += 2) {
                    setpz2(x+p, mulpd2(srN, getpz2(x+p)));
                }
            }
        }
        else
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N; p += 2) {
                setpz2(x+p, mulpd2(srN, getpz2(x+p)));
            }
        }
    }

    void fwdn(complex_vector x) const NOEXCEPT { fwd(x); }

    ///////////////////////////////////////////////////////////////////////////

    void inv(complex_vector x) const NOEXCEPT
    {
        invbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
        }
        else {
            #pragma omp parallel for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
        }
    }

    void inv0(complex_vector x) const NOEXCEPT { inv(x); }

    void invu(complex_vector x) const NOEXCEPT
    {
        const double ssrN = sqrt(1.0/N);
        const ymm srN = { ssrN, ssrN, ssrN, ssrN };
        invbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            if (N >= 2) {
                for (int p = 0; p < N; p += 2)
                    setpz2(x+p, mulpd2(srN, getpz2(x+p)));
            }
        }
        else
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N; p += 2) {
                setpz2(x+p, mulpd2(srN, getpz2(x+p)));
            }
        }
    }

    void invn(complex_vector x) const NOEXCEPT
    {
        const ymm rN = { 1.0/N, 1.0/N, 1.0/N, 1.0/N };
        invbut(N, x, W);
        if (N < OMP_THRESHOLD) {
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            if (N >= 2) {
                for (int p = 0; p < N; p += 2)
                    setpz2(x+p, mulpd2(rN, getpz2(x+p)));
            }
        }
        else
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int p = 0; p < N; p++) {
                const int q = bitrev[p];
                if (p > q) std::swap(x[p], x[q]);
            }
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N; p += 2) {
                setpz2(x+p, mulpd2(rN, getpz2(x+p)));
            }
        }
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // simple_fft2_h
