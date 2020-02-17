/******************************************************************************
*  OTFFT C & Fortran Interface Version 11.5e
*
*  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_c_cpp
#define otfft_c_cpp

extern "C" {
    //=========================================================================

    void *simd_malloc(size_t n) { return OTFFT::simd_malloc(n); }
    void simd_free(void *p) { OTFFT::simd_free(p); }

    //=========================================================================

    void *otfft_fft_new(int N)
    {
        OTFFT::FFT *fft = 0;
        try { fft = new OTFFT::FFT(N); }
        catch (...) { return 0; }
        return fft;
    }
    void otfft_fft_delete(void *p)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        delete fft;
    }

    void otfft_fft_fwd(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->fwd(x_);
    }
    void otfft_fft_fwd0(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->fwd0(x_);
    }
    void otfft_fft_fwdu(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->fwdu(x_);
    }
    void otfft_fft_fwdn(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->fwdn(x_);
    }

    void otfft_fft_inv(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->inv(x_);
    }
    void otfft_fft_inv0(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->inv0(x_);
    }
    void otfft_fft_invu(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->invu(x_);
    }
    void otfft_fft_invn(void *p, double _Complex *x)
    {
        OTFFT::FFT *fft = (OTFFT::FFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        fft->invn(x_);
    }

    //=========================================================================

    void *otfft_fft0_new(int N)
    {
        OTFFT::FFT0 *fft = 0;
        try { fft = new OTFFT::FFT0(N); }
        catch (...) { return 0; }
        return fft;
    }
    void otfft_fft0_delete(void *p)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        delete fft;
    }

    void otfft_fft0_fwd(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->fwd(x_, y_);
    }
    void otfft_fft0_fwd0(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->fwd0(x_, y_);
    }
    void otfft_fft0_fwdu(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->fwdu(x_, y_);
    }
    void otfft_fft0_fwdn(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->fwdn(x_, y_);
    }

    void otfft_fft0_inv(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->inv(x_, y_);
    }
    void otfft_fft0_inv0(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->inv0(x_, y_);
    }
    void otfft_fft0_invu(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->invu(x_, y_);
    }
    void otfft_fft0_invn(void *p, double _Complex *x, double _Complex *y)
    {
        OTFFT::FFT0 *fft = (OTFFT::FFT0 *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        fft->invn(x_, y_);
    }

    //=========================================================================

    void *otfft_rfft_new(int N)
    {
        OTFFT::RFFT *rfft = 0;
        try { rfft = new OTFFT::RFFT(N); }
        catch (...) { return 0; }
        return rfft;
    }
    void otfft_rfft_delete(void *p)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        delete rfft;
    }

    void otfft_rfft_fwd(void *p, const double *x, double _Complex *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        rfft->fwd(x, y_);
    }
    void otfft_rfft_fwd0(void *p, const double *x, double _Complex *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        rfft->fwd0(x, y_);
    }
    void otfft_rfft_fwdu(void *p, const double *x, double _Complex *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        rfft->fwdu(x, y_);
    }
    void otfft_rfft_fwdn(void *p, const double *x, double _Complex *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *y_ = (OTFFT::complex_t *) y;
        rfft->fwdn(x, y_);
    }

    void otfft_rfft_inv(void *p, double _Complex *x, double *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        rfft->inv(x_, y);
    }
    void otfft_rfft_inv0(void *p, double _Complex *x, double *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        rfft->inv0(x_, y);
    }
    void otfft_rfft_invu(void *p, double _Complex *x, double *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        rfft->invu(x_, y);
    }
    void otfft_rfft_invn(void *p, double _Complex *x, double *y)
    {
        OTFFT::RFFT *rfft = (OTFFT::RFFT *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        rfft->invn(x_, y);
    }

    //=========================================================================

    void *otfft_dct_new(int N)
    {
        OTFFT::DCT *dct = 0;
        try { dct = new OTFFT::DCT(N); }
        catch (...) { return 0; }
        return dct;
    }
    void otfft_dct_delete(void *p)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        delete dct;
    }

    void otfft_dct_fwd(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->fwd(x);
    }
    void otfft_dct_fwd0(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->fwd0(x);
    }
    void otfft_dct_fwdn(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->fwdn(x);
    }

    void otfft_dct_inv(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->inv(x);
    }
    void otfft_dct_inv0(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->inv0(x);
    }
    void otfft_dct_invn(void *p, double *x)
    {
        OTFFT::DCT *dct = (OTFFT::DCT *) p;
        dct->invn(x);
    }

    //=========================================================================

    void *otfft_dct0_new(int N)
    {
        OTFFT::DCT0 *dct = 0;
        try { dct = new OTFFT::DCT0(N); }
        catch (...) { return 0; }
        return dct;
    }
    void otfft_dct0_delete(void *p)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        delete dct;
    }

    void otfft_dct0_fwd(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->fwd(x, y, z_);
    }
    void otfft_dct0_fwd0(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->fwd0(x, y, z_);
    }
    void otfft_dct0_fwdn(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->fwdn(x, y, z_);
    }

    void otfft_dct0_inv(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->inv(x, y, z_);
    }
    void otfft_dct0_inv0(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->inv0(x, y, z_);
    }
    void otfft_dct0_invn(void *p, double *x, double *y, double _Complex *z)
    {
        OTFFT::DCT0 *dct = (OTFFT::DCT0 *) p;
        OTFFT::complex_t *z_ = (OTFFT::complex_t *) z;
        dct->invn(x, y, z_);
    }

    //=========================================================================

    void *otfft_bluestein_new(int N)
    {
        OTFFT::Bluestein *bst = 0;
        try { bst = new OTFFT::Bluestein(N); }
        catch (...) { return 0; }
        return bst;
    }
    void otfft_bluestein_delete(void *p)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        delete bst;
    }

    void otfft_bluestein_fwd(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->fwd(x_);
    }
    void otfft_bluestein_fwd0(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->fwd0(x_);
    }
    void otfft_bluestein_fwdu(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->fwdu(x_);
    }
    void otfft_bluestein_fwdn(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->fwdn(x_);
    }

    void otfft_bluestein_inv(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->inv(x_);
    }
    void otfft_bluestein_inv0(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->inv0(x_);
    }
    void otfft_bluestein_invu(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->invu(x_);
    }
    void otfft_bluestein_invn(void *p, double _Complex *x)
    {
        OTFFT::Bluestein *bst = (OTFFT::Bluestein *) p;
        OTFFT::complex_t *x_ = (OTFFT::complex_t *) x;
        bst->invn(x_);
    }

    //=========================================================================
}

#endif // otfft_c.cpp
