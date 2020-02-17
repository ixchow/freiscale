/******************************************************************************
*  OTFFT Header Version 11.5e
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_h
#define otfft_h

#ifdef __cplusplus

#include "otfft_complex.h"

namespace OTFFT { /////////////////////////////////////////////////////////////

using namespace OTFFT_Complex;

/******************************************************************************
*  Complex FFT
******************************************************************************/

struct FFT0
{
    void* obj;
    int N, log_N;

    FFT0() NOEXCEPT;
    FFT0(int n);
    ~FFT0() NOEXCEPT;

    void setup(int n);

    void fwd(complex_vector  x, complex_vector y) const NOEXCEPT;
    void fwd0(complex_vector x, complex_vector y) const NOEXCEPT;
    void fwdu(complex_vector x, complex_vector y) const NOEXCEPT;
    void fwdn(complex_vector x, complex_vector y) const NOEXCEPT;
    void inv(complex_vector  x, complex_vector y) const NOEXCEPT;
    void inv0(complex_vector x, complex_vector y) const NOEXCEPT;
    void invu(complex_vector x, complex_vector y) const NOEXCEPT;
    void invn(complex_vector x, complex_vector y) const NOEXCEPT;
};

struct FFT
{
    FFT0 fft;
    simd_array<complex_t> work;
    complex_t* y;

    FFT() NOEXCEPT : fft(), work(), y(0) {}
    FFT(int n) : fft(n), work(n), y(&work) {}

    inline void setup(int n) { fft.setup(n); work.setup(n); y = &work; }

    inline void fwd(complex_vector  x) const NOEXCEPT { fft.fwd(x, y);  }
    inline void fwd0(complex_vector x) const NOEXCEPT { fft.fwd0(x, y); }
    inline void fwdu(complex_vector x) const NOEXCEPT { fft.fwdu(x, y); }
    inline void fwdn(complex_vector x) const NOEXCEPT { fft.fwdn(x, y); }
    inline void inv(complex_vector  x) const NOEXCEPT { fft.inv(x, y);  }
    inline void inv0(complex_vector x) const NOEXCEPT { fft.inv0(x, y); }
    inline void invu(complex_vector x) const NOEXCEPT { fft.invu(x, y); }
    inline void invn(complex_vector x) const NOEXCEPT { fft.invn(x, y); }
};

/******************************************************************************
*  Real FFT
******************************************************************************/

struct RFFT
{
#ifdef DO_SINGLE_THREAD
    static const int OMP_THRESHOLD   = 1<<30;
    static const int OMP_THRESHOLD_W = 1<<30;
#else
    static const int OMP_THRESHOLD   = 1<<15;
    static const int OMP_THRESHOLD_W = 1<<16;
#endif

    int N;
    FFT0 fft;
    simd_array<complex_t> weight;
    complex_t* U;

    RFFT() NOEXCEPT;
    RFFT(int n);

    void setup(int n);

    void fwd(const_double_vector  x, complex_vector y) const NOEXCEPT;
    void fwd0(const_double_vector x, complex_vector y) const NOEXCEPT;
    void fwdu(const_double_vector x, complex_vector y) const NOEXCEPT;
    void fwdn(const_double_vector x, complex_vector y) const NOEXCEPT;
    void inv(complex_vector  x, double_vector y) const NOEXCEPT;
    void inv0(complex_vector x, double_vector y) const NOEXCEPT;
    void invu(complex_vector x, double_vector y) const NOEXCEPT;
    void invn(complex_vector x, double_vector y) const NOEXCEPT;
};

/******************************************************************************
*  DCT
******************************************************************************/

struct DCT0
{
#ifdef DO_SINGLE_THREAD
    static const int OMP_THRESHOLD   = 1<<30;
    static const int OMP_THRESHOLD_W = 1<<30;
#else
    static const int OMP_THRESHOLD   = 1<<15;
    static const int OMP_THRESHOLD_W = 1<<16;
#endif

    int N;
    RFFT rfft;
    simd_array<complex_t> weight;
    complex_t* V;

    DCT0() NOEXCEPT;
    DCT0(int n);

    void setup(int n);

    void fwd(double_vector  x, double_vector y, complex_vector z) const NOEXCEPT;
    void fwd0(double_vector x, double_vector y, complex_vector z) const NOEXCEPT;
    void fwdn(double_vector x, double_vector y, complex_vector z) const NOEXCEPT;
    void inv(double_vector  x, double_vector y, complex_vector z) const NOEXCEPT;
    void inv0(double_vector x, double_vector y, complex_vector z) const NOEXCEPT;
    void invn(double_vector x, double_vector y, complex_vector z) const NOEXCEPT;
};

struct DCT
{
    int N;
    DCT0 dct;
    simd_array<double> work1;
    simd_array<complex_t> work2;
    double* y;
    complex_t* z;

    DCT() NOEXCEPT : N(0), y(0), z(0) {}
    DCT(int n) { setup(n); }

    inline void setup(int n)
    {
        N = n;
        dct.setup(N);
        work1.setup(N); y = &work1;
        work2.setup(N); z = &work2;
    }

    inline void fwd(double_vector  x) const NOEXCEPT { dct.fwd(x, y, z);  }
    inline void fwd0(double_vector x) const NOEXCEPT { dct.fwd0(x, y, z); }
    inline void fwdn(double_vector x) const NOEXCEPT { dct.fwdn(x, y, z); }
    inline void inv(double_vector  x) const NOEXCEPT { dct.inv(x, y, z);  }
    inline void inv0(double_vector x) const NOEXCEPT { dct.inv0(x, y, z); }
    inline void invn(double_vector x) const NOEXCEPT { dct.invn(x, y, z); }
};

/******************************************************************************
*  Bluestein's FFT
******************************************************************************/

struct Bluestein
{
#ifdef DO_SINGLE_THREAD
    static const int OMP_THRESHOLD   = 1<<30;
    static const int OMP_THRESHOLD_W = 1<<30;
#else
    static const int OMP_THRESHOLD   = 1<<15;
    static const int OMP_THRESHOLD_W = 1<<16;
#endif

    int N, L;
    FFT fft;
    simd_array<complex_t> work1;
    simd_array<complex_t> work2;
    simd_array<complex_t> weight;
    complex_t* a;
    complex_t* b;
    complex_t* W;

    Bluestein() NOEXCEPT;
    Bluestein(int n);

    void setup(int n);

    void fwd(complex_vector  x) const NOEXCEPT;
    void fwd0(complex_vector x) const NOEXCEPT;
    void fwdu(complex_vector x) const NOEXCEPT;
    void fwdn(complex_vector x) const NOEXCEPT;
    void inv(complex_vector  x) const NOEXCEPT;
    void inv0(complex_vector x) const NOEXCEPT;
    void invu(complex_vector x) const NOEXCEPT;
    void invn(complex_vector x) const NOEXCEPT;
};

} /////////////////////////////////////////////////////////////////////////////

#else
#include "otfft_c.h"
#endif // __cplusplus

#endif // otfft_h
