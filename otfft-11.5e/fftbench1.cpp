/******************************************************************************
*  FFT Benchmark 1
******************************************************************************/

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include "cpp_fftw3.h"
#include "ooura1.h"
#include "simple_fft.h"
//#include "simple_fft2.h"
#include "otfft/otfft.h"
#include "otfft/msleep.h"

typedef std::chrono::microseconds::rep counter_t;

using namespace std;
using OTFFT::complex_t;
using OTFFT::simd_malloc;
using OTFFT::simd_free;
using CppFFTW3::FORWARD;
using CppFFTW3::INVERSE;

#define DELAY1 1
#define DELAY2 100
#define FACTOR 1

double safe_avr(const vector<counter_t>& dt)
{
    typedef vector<counter_t>::size_type size_type;
    const size_type TRIES = dt.size();
    counter_t sum = 0;
    for (size_type i = 0; i < TRIES; i++) sum += dt[i];
    const double m = double(sum)/TRIES;
    double sum_dd = 0;
    for (size_type i = 0; i < TRIES; i++) {
        const double d = dt[i] - m;
        sum_dd += d*d;
    }
    const double ss = sum_dd/TRIES;
    sum = 0;
    size_type n = 0;
    for (size_type i = 0; i < TRIES; i++) {
        const double d = dt[i] - m;
        if (d*d <= 2*ss) { sum += dt[i]; n++; }
        //if (d*d <= 1*ss) { sum += dt[i]; n++; }
    }
    return double(sum)/n;
}

template <typename FFT, typename IFFT>
double laptime1(int LOOPS, int TRIES, const FFT& fft, const IFFT& ifft)
{
    using namespace chrono;
    vector<counter_t> dt(TRIES);
    for (int i = 0; i < TRIES; i++) {
        const system_clock::time_point t1 = system_clock::now();
        for (int j = 0; j < LOOPS; j++) {
            fft();
            ifft();
        }
        const system_clock::time_point t2 = system_clock::now();
        dt[i] = duration_cast<microseconds>(t2 - t1).count();
        msleep(DELAY1);
    }
    return safe_avr(dt);
}

template <typename FFT>
double laptime2(int LOOPS, int TRIES, const FFT& fft, complex_t *x)
{
    using namespace chrono;
    vector<counter_t> dt(TRIES);
    for (int i = 0; i < TRIES; i++) {
        const system_clock::time_point t1 = system_clock::now();
        for (int j = 0; j < LOOPS; j++) {
            fft.fwd(x);
            fft.inv(x);
        }
        const system_clock::time_point t2 = system_clock::now();
        dt[i] = duration_cast<microseconds>(t2 - t1).count();
        msleep(DELAY1);
    }
    return safe_avr(dt);
}

struct fd {
    const int n;
    const int x;
    fd(const int n, const int x) : n(n), x(x) {}
};
ostream& operator<<(ostream& os, const fd& f)
{
    os << setw(f.n) << f.x;
    return os;
}

struct ff {
    const int n;
    const int m;
    const double x;
    ff(const int n, const int m, const double x) : n(n), m(m), x(x) {}
};
ostream& operator<<(ostream& os, const ff& f)
{
    os << fixed << setw(f.n) << setprecision(f.m) << f.x;
    return os;
}

void initialize(const int N, complex_t* x)
{
    for (int p = 0; p < N; p++) {
        const double t = double(p)/N;
        x[p].Re = 10 * cos(3*2*M_PI*t*t);
        x[p].Im = 10 * sin(3*2*M_PI*t*t);
    }
}

int main(int argc, char *argv[]) try
{
    const int n_min  = argc >= 2 ? atoi(argv[1]) : 1;
    const int n_max  = argc >= 3 ? atoi(argv[2]) : (argc == 2 ? n_min : 22);
    const int N_max  = 1 << n_max;
    const int nN_max = (max)(22*(1<<22), n_max*N_max);
    if (n_min < 1 || 24 < n_min) throw "argv[1] must be 1..24";
    if (n_max < 1 || 24 < n_max) throw "argv[2] must be 1..24";
    if (n_min > n_max) throw "argv must be argv[1] <= argv[2]";

#ifndef DO_SINGLE_THREAD
    fftw_init_threads();
#endif
    cout << "------+-----------+-----------------+-----------------+-----------------+---\n";
    cout << "length|FFTW3[usec]|   OOURA   [usec]| SimpleFFT [usec]|   OTFFT   [usec]|err\n";
    cout << "------+-----------+-----------------+-----------------+-----------------+---\n";
    complex_t* x1 = (complex_t*) simd_malloc(N_max*sizeof(complex_t));
    complex_t* x2 = (complex_t*) simd_malloc(N_max*sizeof(complex_t));
    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    for (int n = n_min; n <= n_max; n++) {
        const int N = 1 << n;
        const int LOOPS = nN_max/(n*N);
        const int TRIES = (min)(16, n*FACTOR);
        //const int TRIES = (min)(64, n*FACTOR);
        double lap, lap1;

        cout << "2^(" << fd(2, n) << ")|" << flush;

        //// FFTW3 ////
        initialize(N, x1);
        CppFFTW3::Plan<FORWARD,1,1> fwd_plan(N, x1);
        CppFFTW3::Plan<INVERSE,0,0> inv_plan(N, x1);
        lap = laptime1(LOOPS, TRIES, fwd_plan, inv_plan);
        sum1 += lap*n/n_max;
        cout << ff(11,2, lap/LOOPS) << "|" << flush;
        msleep(DELAY2);
        lap1 = lap;

        //// OOURA FFT ////
        initialize(N, x1);
        lap = laptime2(LOOPS, TRIES, OOURA::FFT(N), x1);
        sum2 += lap*n/n_max;
        cout << ff(11,2, lap/LOOPS) << "(" << ff(3,0, 100*lap/lap1) << "%)|" << flush;
        msleep(DELAY2);

        // erase OOURA FFT's influence
        CppFFTW3::Plan<FORWARD,1,1> eraser(N, x1); eraser();

        //// Simple FFT ////
        initialize(N, x1);
        //OTFFT::speedup_magic();
        lap = laptime2(LOOPS, TRIES, SimpleFFT::FFT(N), x1);
        sum3 += lap*n/n_max;
        cout << ff(11,2, lap/LOOPS) << "(" << ff(3,0, 100*lap/lap1) << "%)|" << flush;
        msleep(DELAY2);

        //// OTFFT ////
        initialize(N, x1);
        //OTFFT::speedup_magic();
        lap = laptime2(LOOPS, TRIES, OTFFT::FFT(N), x1);
        sum4 += lap*n/n_max;
        cout << ff(11,2, lap/LOOPS) << "(" << ff(3,0, 100*lap/lap1) << "%)|" << flush;
        msleep(DELAY2);

        double err = 0;
        initialize(N, x1);
        initialize(N, x2);
        fwd_plan();
        //SimpleFFT::FFT(N).fwd(x1);
        OTFFT::FFT(N).fwd(x2);
        for (int k = 0; k < N; k++) {
            const complex_t d = x2[k] - x1[k];
            err += Re(d*conj(d));
        }
        initialize(N, x1);
        initialize(N, x2);
        inv_plan();
        //SimpleFFT::FFT(N).inv(x1);
        OTFFT::FFT(N).inv(x2);
        for (int k = 0; k < N; k++) {
            const complex_t d = x2[k] - x1[k];
            err += Re(d*conj(d));
        }
        if (err == 0)
            cout << " -" << endl;
        else if (err > 0)
            cout << ff(3,0, log10(err)) << endl;
        else
            cout << "NG" << endl;
        msleep(DELAY2);
    }
    simd_free(x2);
    simd_free(x1);
    cout << "------+-----------+-----------------+-----------------+-----------------+---\n";
    cout << " cost |";
    cout << ff(11,2, sum1/(n_max-n_min+1)) << "|" << flush;
    cout << ff(11,2, sum2/(n_max-n_min+1)) << "(" << ff(3,0, 100*sum2/sum1) << "%)|" << flush;
    cout << ff(11,2, sum3/(n_max-n_min+1)) << "(" << ff(3,0, 100*sum3/sum1) << "%)|" << flush;
    cout << ff(11,2, sum4/(n_max-n_min+1)) << "(" << ff(3,0, 100*sum4/sum1) << "%)|" << endl;
    cout << "------+-----------+-----------------+-----------------+-----------------+---\n";
#ifndef DO_SINGLE_THREAD
    fftw_cleanup_threads();
#endif

    return 0;
}
catch (const char *message) { cerr << message << endl; exit(1); }
catch (...) { cerr << "\n""*** exception! ***" << endl; exit(1); }
