/******************************************************************************
*  FFT Tuning Command Version 11.5e
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include "otfft_misc.h"
#include "otfft_avxdif4.h"
#include "otfft_avxdit4.h"
#include "otfft_avxdif8.h"
#include "otfft_avxdit8.h"
#include "otfft_avxdif16.h"
#include "otfft_avxdit16.h"
#include "otfft_sixstep.h"
#include "msleep.h"

#ifdef _MSC_VER
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;
using namespace OTFFT_MISC;

#define DELAY1 1
#define DELAY2 100
#define FACTOR 1

static const int n_max = 24;

template <class FFT>
double laptime(int LOOPS, int TRIES, const FFT& fft, complex_t *x, complex_t *y)
{
    using namespace chrono;
    typedef microseconds::rep counter_t;
    static const int th = (n_max*(1<<n_max))/(14*(1<<14));
    counter_t sum = 0;
    vector<counter_t> dt(TRIES);
    for (int i = 0; i < TRIES; i++) {
        //speedup_magic();
        const system_clock::time_point t1 = system_clock::now();
        for (int j = 0; j < LOOPS; j++) {
            fft.fwd(x, y);
            fft.inv(x, y);
        }
        const system_clock::time_point t2 = system_clock::now();
        sum += (dt[i] = duration_cast<microseconds>(t2 - t1).count());
        msleep(DELAY1);
    }
    const double m = double(sum)/TRIES;
    double sum_dd = 0;
    for (int i = 0; i < TRIES; i++) {
        const double d = dt[i] - m;
        sum_dd += d*d;
    }
    const double ss = sum_dd/TRIES;
    sum = 0;
    int n = 0;
    for (int i = 0; i < TRIES; i++) {
        const double d = dt[i] - m;
        if (d*d <= 2*ss) { sum += dt[i]; n++; }
    }
    if (LOOPS >= th)
        cout << fixed << setw(10) << setprecision(3) << double(sum)/n/LOOPS << flush;
    else
        cout << fixed << setw(10) << setprecision(0) << double(sum)/n/LOOPS << flush;
    return double(sum)/n;
}

void initialize(const int N, complex_vector x)
{
    for (int p = 0; p < N; p++) {
        const double t = double(p)/N;
        x[p].Re = 10 * cos(3*2*M_PI*t*t);
        x[p].Im = 10 * sin(3*2*M_PI*t*t);
    }
}

int tune(int n)
{
    static const int N_max  = 1 << n_max;
    static const int nN_max = n_max * N_max;
    const int N = 1 << n;
    const int LOOPS = nN_max/(n*N);
    const int TRIES = (min)(16, n*FACTOR);
    double lap, tmp;
    int fastest = 1;
    complex_vector x = (complex_vector) simd_malloc(N*sizeof(complex_t));
    complex_vector y = (complex_vector) simd_malloc(N*sizeof(complex_t));

    cout << "2^(" << setw(2) << n << ")" << flush;

    initialize(N, x);
    lap = laptime(LOOPS, TRIES, OTFFT_AVXDIF4::FFT0(N), x, y);

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_AVXDIT4::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 2; }

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_AVXDIF8::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 3; }

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_AVXDIT8::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 4; }

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_AVXDIF16::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 5; }

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_AVXDIT16::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 6; }

    msleep(DELAY2);

    initialize(N, x);
    tmp = laptime(LOOPS, TRIES, OTFFT_SixStep::FFT0(N), x, y);
    if (tmp < lap) { lap = tmp; fastest = 7; }

    msleep(DELAY2);
    cout << " :" << fastest << endl;

    simd_free(y);
    simd_free(x);

    return fastest;
}

void generate()
{
    ofstream ofs0("otfft_gen_new.h");
    ofstream ofs1("otfft_gen_setup.h");
    ofstream ofs2("otfft_gen_fwd.h");
    ofstream ofs3("otfft_gen_inv.h");
    ofstream ofs4("otfft_gen_fwd0.h");
    ofstream ofs5("otfft_gen_invn.h");
    ofstream ofs6("otfft_gen_fwdu.h");
    ofstream ofs7("otfft_gen_invu.h");
    ofstream ofs8("otfft_gen_delete.h");
    ofs0 << "switch (log_N) {\ncase  0: break;\n";
    ofs1 << "switch (log_N) {\ncase  0: break;\n";
    ofs2 << "switch (log_N) {\ncase  0: break;\n";
    ofs3 << "switch (log_N) {\ncase  0: break;\n";
    ofs4 << "switch (log_N) {\ncase  0: break;\n";
    ofs5 << "switch (log_N) {\ncase  0: break;\n";
    ofs6 << "switch (log_N) {\ncase  0: break;\n";
    ofs7 << "switch (log_N) {\ncase  0: break;\n";
    ofs8 << "switch (log_N) {\ncase  0: break;\n";
    for (int n = 1; n <= n_max; n++) {
        const int fastest = tune(n);

        ofs0 << "case " << setw(2) << n << ": obj = new FFT" << fastest << "(); break;\n";
        ofs1 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->setup2(log_N); break;\n";
        ofs2 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwd(x, y); break;\n";
        ofs3 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->inv(x, y); break;\n";
        ofs4 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwd0(x, y); break;\n";
        ofs5 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->invn(x, y); break;\n";
        ofs6 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwdu(x, y); break;\n";
        ofs7 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->invu(x, y); break;\n";
        ofs8 << "case " << setw(2) << n << ": delete static_cast<FFT" << fastest << "*>(obj); break;\n";
    }
    ofs0 << "default: obj = new FFT8(); break;\n";
    ofs1 << "default: static_cast<FFT8*>(obj)->setup(N); break;\n";
    ofs2 << "default: static_cast<FFT8*>(obj)->fwd(x, y); break;\n";
    ofs3 << "default: static_cast<FFT8*>(obj)->inv(x, y); break;\n";
    ofs4 << "default: static_cast<FFT8*>(obj)->fwd0(x, y); break;\n";
    ofs5 << "default: static_cast<FFT8*>(obj)->invn(x, y); break;\n";
    ofs6 << "default: static_cast<FFT8*>(obj)->fwdu(x, y); break;\n";
    ofs7 << "default: static_cast<FFT8*>(obj)->invu(x, y); break;\n";
    ofs8 << "default: delete static_cast<FFT8*>(obj); break;\n";
    ofs0 << "}\n";
    ofs1 << "}\n";
    ofs2 << "}\n";
    ofs3 << "}\n";
    ofs4 << "}\n";
    ofs5 << "}\n";
    ofs6 << "}\n";
    ofs7 << "}\n";
    ofs8 << "}\n";
}

void oneline(const int n)
{
    string line[9][n_max];
    {
        ifstream ifs0("otfft_gen_new.h");
        ifstream ifs1("otfft_gen_setup.h");
        ifstream ifs2("otfft_gen_fwd.h");
        ifstream ifs3("otfft_gen_inv.h");
        ifstream ifs4("otfft_gen_fwd0.h");
        ifstream ifs5("otfft_gen_invn.h");
        ifstream ifs6("otfft_gen_fwdu.h");
        ifstream ifs7("otfft_gen_invu.h");
        ifstream ifs8("otfft_gen_delete.h");
        string header;
        getline(ifs0, header); getline(ifs0, header);
        getline(ifs1, header); getline(ifs1, header);
        getline(ifs2, header); getline(ifs2, header);
        getline(ifs3, header); getline(ifs3, header);
        getline(ifs4, header); getline(ifs4, header);
        getline(ifs5, header); getline(ifs5, header);
        getline(ifs6, header); getline(ifs6, header);
        getline(ifs7, header); getline(ifs7, header);
        getline(ifs8, header); getline(ifs8, header);
        for (int i = 0; i < n_max; i++) getline(ifs0, line[0][i]);
        for (int i = 0; i < n_max; i++) getline(ifs1, line[1][i]);
        for (int i = 0; i < n_max; i++) getline(ifs2, line[2][i]);
        for (int i = 0; i < n_max; i++) getline(ifs3, line[3][i]);
        for (int i = 0; i < n_max; i++) getline(ifs4, line[4][i]);
        for (int i = 0; i < n_max; i++) getline(ifs5, line[5][i]);
        for (int i = 0; i < n_max; i++) getline(ifs6, line[6][i]);
        for (int i = 0; i < n_max; i++) getline(ifs7, line[7][i]);
        for (int i = 0; i < n_max; i++) getline(ifs8, line[8][i]);
    }

    ofstream ofs0("otfft_gen_new.h");
    ofstream ofs1("otfft_gen_setup.h");
    ofstream ofs2("otfft_gen_fwd.h");
    ofstream ofs3("otfft_gen_inv.h");
    ofstream ofs4("otfft_gen_fwd0.h");
    ofstream ofs5("otfft_gen_invn.h");
    ofstream ofs6("otfft_gen_fwdu.h");
    ofstream ofs7("otfft_gen_invu.h");
    ofstream ofs8("otfft_gen_delete.h");
    ofs0 << "switch (log_N) {\ncase  0: break;\n";
    ofs1 << "switch (log_N) {\ncase  0: break;\n";
    ofs2 << "switch (log_N) {\ncase  0: break;\n";
    ofs3 << "switch (log_N) {\ncase  0: break;\n";
    ofs4 << "switch (log_N) {\ncase  0: break;\n";
    ofs5 << "switch (log_N) {\ncase  0: break;\n";
    ofs6 << "switch (log_N) {\ncase  0: break;\n";
    ofs7 << "switch (log_N) {\ncase  0: break;\n";
    ofs8 << "switch (log_N) {\ncase  0: break;\n";

    const int fastest = tune(n);
    ostringstream oss0, oss1, oss2, oss3, oss4, oss5, oss6, oss7, oss8;
    oss0 << "case " << setw(2) << n << ": obj = new FFT" << fastest << "(); break;";
    oss1 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->setup2(log_N); break;";
    oss2 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwd(x, y); break;";
    oss3 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->inv(x, y); break;";
    oss4 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwd0(x, y); break;";
    oss5 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->invn(x, y); break;";
    oss6 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->fwdu(x, y); break;";
    oss7 << "case " << setw(2) << n << ": static_cast<FFT" << fastest << "*>(obj)->invu(x, y); break;";
    oss8 << "case " << setw(2) << n << ": delete static_cast<FFT" << fastest << "*>(obj); break;";
    line[0][n-1] = oss0.str();
    line[1][n-1] = oss1.str();
    line[2][n-1] = oss2.str();
    line[3][n-1] = oss3.str();
    line[4][n-1] = oss4.str();
    line[5][n-1] = oss5.str();
    line[6][n-1] = oss6.str();
    line[7][n-1] = oss7.str();
    line[8][n-1] = oss8.str();
    for (int i = 0; i < n_max; i++) ofs0 << line[0][i] << endl;
    for (int i = 0; i < n_max; i++) ofs1 << line[1][i] << endl;
    for (int i = 0; i < n_max; i++) ofs2 << line[2][i] << endl;
    for (int i = 0; i < n_max; i++) ofs3 << line[3][i] << endl;
    for (int i = 0; i < n_max; i++) ofs4 << line[4][i] << endl;
    for (int i = 0; i < n_max; i++) ofs5 << line[5][i] << endl;
    for (int i = 0; i < n_max; i++) ofs6 << line[6][i] << endl;
    for (int i = 0; i < n_max; i++) ofs7 << line[7][i] << endl;
    for (int i = 0; i < n_max; i++) ofs8 << line[8][i] << endl;

    ofs0 << "default: obj = new FFT8(); break;\n";
    ofs1 << "default: static_cast<FFT8*>(obj)->setup(N); break;\n";
    ofs2 << "default: static_cast<FFT8*>(obj)->fwd(x, y); break;\n";
    ofs3 << "default: static_cast<FFT8*>(obj)->inv(x, y); break;\n";
    ofs4 << "default: static_cast<FFT8*>(obj)->fwd0(x, y); break;\n";
    ofs5 << "default: static_cast<FFT8*>(obj)->invn(x, y); break;\n";
    ofs6 << "default: static_cast<FFT8*>(obj)->fwdu(x, y); break;\n";
    ofs7 << "default: static_cast<FFT8*>(obj)->invu(x, y); break;\n";
    ofs8 << "default: delete static_cast<FFT8*>(obj); break;\n";
    ofs0 << "}\n";
    ofs1 << "}\n";
    ofs2 << "}\n";
    ofs3 << "}\n";
    ofs4 << "}\n";
    ofs5 << "}\n";
    ofs6 << "}\n";
    ofs7 << "}\n";
    ofs8 << "}\n";
}

void cd(const char* dir)
{
#if _MSC_VER
    const int err = _chdir(dir);
#else
    const int err = chdir(dir);
#endif
    if (err != 0) throw "chdir: failed";
}

int main(int argc, char *argv[]) try
{
    if (argc == 1) {
        cout << "Optimizing from 2^( 1) to 2^(" << n_max << ")" << endl;
        generate();
        cout << "Done!" << endl;
    }
    else if (argc == 2) {
        const int n = atoi(argv[1]);
        if (1 <= n && n <= n_max) oneline(n);
        else throw "length: out of range";
    }
    else if (argc == 3 && argv[1][0] == '@') {
        cd(argv[2]);
        cout << "Optimizing from 2^( 1) to 2^(" << n_max << ")" << endl;
        generate();
        cout << "Done!" << endl;
    }
    else if (argc == 4 && argv[1][0] == '@') {
        cd(argv[2]);
        const int n = atoi(argv[3]);
        if (1 <= n && n <= n_max) oneline(n);
        else throw "length: out of range";
    }
    else throw "usage1: ffttune [length]\n"
               "usage2: ffttune @ otfft_path [length]";
    return 0;
}
catch (const char* message) { cerr << message << endl; exit(1); }
catch (...) { cerr << "\n""*** exception! ***\n"; exit(1); }
