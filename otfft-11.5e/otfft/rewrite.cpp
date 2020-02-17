/******************************************************************************
*  Rewrite Command Version 11.5e
*
*  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

constexpr int n_max = 24;

void rewrite(const int n, const int fft_num)
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

    ostringstream oss0, oss1, oss2, oss3, oss4, oss5, oss6, oss7, oss8;
    oss0 << "case " << setw(2) << n << ": obj = new FFT" << fft_num << "(); break;";
    oss1 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->setup2(log_N); break;";
    oss2 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->fwd(x, y); break;";
    oss3 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->inv(x, y); break;";
    oss4 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->fwd0(x, y); break;";
    oss5 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->invn(x, y); break;";
    oss6 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->fwdu(x, y); break;";
    oss7 << "case " << setw(2) << n << ": static_cast<FFT" << fft_num << "*>(obj)->invu(x, y); break;";
    oss8 << "case " << setw(2) << n << ": delete static_cast<FFT" << fft_num << "*>(obj); break;";
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

int main(int argc, char *argv[]) try
{
    const int n = argc == 3 ? atoi(argv[1]) : 0;
    const int m = argc == 3 ? atoi(argv[2]) : 0;
    if (1 <= n && n <= n_max) {
        if (1 <= m && m <= 8) rewrite(n, m);
    }
    else throw "usage: rewrite length fft-number";
    return 0;
}
catch (const char* message) { cerr << message << endl; exit(1); }
catch (...) { cerr << "\n""*** exception! ***\n"; exit(1); }
