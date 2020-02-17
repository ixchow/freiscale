/******************************************************************************
*  OTFFT SixStep Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_sixstep_h
#define otfft_sixstep_h

#include "otfft_misc.h"

namespace OTFFT_SixStep { /////////////////////////////////////////////////////

using namespace OTFFT_MISC;

typedef const_complex_vector weight_t;
struct index_t { short col, row; };
typedef const index_t* __restrict const const_index_vector;

} /////////////////////////////////////////////////////////////////////////////

#include "otfft_sixstepsq.h"
#include "otfft_eightstep.h"

namespace OTFFT_SixStep { /////////////////////////////////////////////////////

using namespace OTFFT_EightStep;

struct FFT0
{
    typedef complex_t* __restrict complex_ptr;
    typedef index_t* __restrict index_ptr;
    int N, log_N;
    complex_ptr W;
    complex_ptr Wm;
    complex_ptr Ws;
    index_ptr iv;
    simd_array<complex_t> weight;
    simd_array<complex_t> weight_med;
    simd_array<complex_t> weight_sub;
    simd_array<index_t> index_array;

    FFT0() noexcept : N(0), log_N(0), W(0), Wm(0), Ws(0), iv(0) {}
    FFT0(const int n) { setup(n); }

    void setup(int n)
    {
        W = 0; Wm = 0; Ws = 0; iv = 0;
        for (log_N = 0; n > 1; n >>= 1) log_N++;
        setup2(log_N);
    }

    void setup2(const int n)
    {
        log_N = n; N = 1 << n;
        if (n <= 5)
            init_weight_sub(N, W, weight);
        else if ((n & 1) == 1) {
            const int m = 1 << (n/2-1);
            init_weight_8xm(N, W, weight);
            init_weight_nxn(m, Wm, weight_med, iv, index_array);
            init_weight_sub(m, Ws, weight_sub);
        }
        else {
            const int m = 1 << n/2;
            init_weight_nxn(m, W, weight, iv, index_array);
            init_weight_sub(m, Ws, weight_sub);
        }
    }

    static void init_weight_8xm(const int n, complex_ptr& w, simd_array<complex_t>& a)
    {
        a.setup(n/8);
        w = &a;
        init_W8xm(n, w);
    }

    static void init_weight_nxn(const int n,
            complex_ptr& w, simd_array<complex_t>& a, index_ptr& ip, simd_array<index_t>& ia)
    {
        const int m = n/2;
        a.setup(n*(n+2)/2);
        w = &a;
        ia.setup(m*(m+1)/2);
        ip = &ia;
        init_Wnxn(n, w, ip);
    }

    static void init_weight_sub(const int n, complex_ptr& w, simd_array<complex_t>& a)
    {
#if 0
        if (n <= 4) w = 0;
        else {
            a.setup(2*n);
            w = &a;
            init_Wt(4, n, w);
        }
#elif 1
        if (n <= 8) w = 0;
        else {
            a.setup(2*n);
            w = &a;
            init_Wt(8, n, w);
        }
#else
        if (n <= 16) w = 0;
        else {
            a.setup(2*n);
            w = &a;
            init_Wt(16, n, w);
        }
#endif
    }

    static void init_W8xm(const int n, complex_vector w)
    {
        const int m = n / 8;
        const double theta = -2*M_PI/n;
        if (m < OMP_THRESHOLD1) {
            for (int p = 0; p < m; p++) {
                w[p] = expj(theta * p);
            }
        }
        else {
            #pragma omp parallel for schedule(static)
            for (int p = 0; p < m; p++) {
                w[p] = expj(theta * p);
            }
        }
    }

    static void init_Wnxn(const int n, complex_vector w, index_ptr ip)
    {
        const int N = n*n;
        const int m = n/2;
        const double theta = -2*M_PI/N;
        if (N < OMP_THRESHOLD1) {
            for (int p = 0; p < n; p += 2) {
                for (int k = p; k < n; k += 2) {
                    const int u = p/2;
                    const int v = k/2;
                    const int i = u*m + v - u*(u+1)/2;
                    w[4*i+0] = expj(theta * (k+0)*(p+0));
                    w[4*i+1] = expj(theta * (k+1)*(p+0));
                    w[4*i+2] = expj(theta * (k+0)*(p+1));
                    w[4*i+3] = expj(theta * (k+1)*(p+1));
                    ip[i].row = p;
                    ip[i].col = k;
                }
            }
        }
        else {
            #pragma omp parallel for schedule(guided)
            for (int p = 0; p < n; p += 2) {
                for (int k = p; k < n; k += 2) {
                    const int u = p/2;
                    const int v = k/2;
                    const int i = u*m + v - u*(u+1)/2;
                    w[4*i+0] = expj(theta * (k+0)*(p+0));
                    w[4*i+1] = expj(theta * (k+1)*(p+0));
                    w[4*i+2] = expj(theta * (k+0)*(p+1));
                    w[4*i+3] = expj(theta * (k+1)*(p+1));
                    ip[i].row = p;
                    ip[i].col = k;
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    inline void fwd(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_length;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::fwdfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::fwdfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::fwdfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::fwdfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::fwdfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: fwdffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: fwdfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: fwdffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: fwdfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: fwdffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: fwdfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: fwdffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: fwdfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: fwdffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: fwdfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: fwdffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: fwdfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: fwdffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: fwdfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: fwdffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: fwdfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: fwdffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: fwdfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: fwdffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }

    inline void fwd0(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_1;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::fwdfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::fwdfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::fwdfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::fwdfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::fwdfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: fwdffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: fwdfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: fwdffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: fwdfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: fwdffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: fwdfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: fwdffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: fwdfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: fwdffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: fwdfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: fwdffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: fwdfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: fwdffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: fwdfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: fwdffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: fwdfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: fwdffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: fwdfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: fwdffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }

    inline void fwdu(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_unitary;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::fwdfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::fwdfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::fwdfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::fwdfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::fwdfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: fwdffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: fwdfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: fwdffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: fwdfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: fwdffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: fwdfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: fwdffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: fwdfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: fwdffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: fwdfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: fwdffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: fwdfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: fwdffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: fwdfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: fwdffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: fwdfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: fwdffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: fwdfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: fwdffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }

    inline void fwdn(complex_vector x, complex_vector y) const noexcept { fwd(x, y); }

    inline void inv(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_1;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::invfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::invfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::invfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::invfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::invfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: invffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: invfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: invffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: invfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: invffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: invfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: invffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: invfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: invffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: invfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: invffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: invfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: invffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: invfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: invffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: invfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: invffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: invfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: invffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }

    inline void inv0(complex_vector x, complex_vector y) const noexcept { inv(x, y); }

    inline void invu(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_unitary;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::invfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::invfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::invfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::invfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::invfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: invffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: invfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: invffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: invfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: invffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: invfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: invffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: invfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: invffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: invfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: invffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: invfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: invffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: invfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: invffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: invfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: invffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: invfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: invffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }

    inline void invn(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_length;
        switch (log_N) {
            case  0: break;
            case  1: OTFFT_FourStep::invfft<(1<<1),1,0,mode>()(x, y, W); break;
            case  2: OTFFT_FourStep::invfft<(1<<2),1,0,mode>()(x, y, W); break;
            case  3: OTFFT_FourStep::invfft<(1<<3),1,0,mode>()(x, y, W); break;
            case  4: OTFFT_FourStep::invfft<(1<<4),1,0,mode>()(x, y, W); break;
            case  5: OTFFT_FourStep::invfft<(1<<5),1,0,mode>()(x, y, W); break;
            case  6: invffts< 6,mode>()(iv, x, y, W, Ws); break;
            case  7: invfftr< 7,mode>()(iv, x, y, W, Wm, Ws); break;
            case  8: invffts< 8,mode>()(iv, x, y, W, Ws); break;
            case  9: invfftr< 9,mode>()(iv, x, y, W, Wm, Ws); break;
            case 10: invffts<10,mode>()(iv, x, y, W, Ws); break;
            case 11: invfftr<11,mode>()(iv, x, y, W, Wm, Ws); break;
            case 12: invffts<12,mode>()(iv, x, y, W, Ws); break;
            case 13: invfftr<13,mode>()(iv, x, y, W, Wm, Ws); break;
            case 14: invffts<14,mode>()(iv, x, y, W, Ws); break;
            case 15: invfftr<15,mode>()(iv, x, y, W, Wm, Ws); break;
            case 16: invffts<16,mode>()(iv, x, y, W, Ws); break;
            case 17: invfftr<17,mode>()(iv, x, y, W, Wm, Ws); break;
            case 18: invffts<18,mode>()(iv, x, y, W, Ws); break;
            case 19: invfftr<19,mode>()(iv, x, y, W, Wm, Ws); break;
            case 20: invffts<20,mode>()(iv, x, y, W, Ws); break;
            case 21: invfftr<21,mode>()(iv, x, y, W, Wm, Ws); break;
            case 22: invffts<22,mode>()(iv, x, y, W, Ws); break;
            case 23: invfftr<23,mode>()(iv, x, y, W, Wm, Ws); break;
            case 24: invffts<24,mode>()(iv, x, y, W, Ws); break;
        }
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_sixstep_h
