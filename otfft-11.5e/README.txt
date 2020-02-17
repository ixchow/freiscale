This software is released under the MIT License, see LICENSE.txt.

-------------------------------------------------------------------------------

[What is the OTFFT]

OTFFT is a high-speed FFT library using the Stockham's algorithm and AVX.
In addition, C++ template metaprogramming technique is used in OTFFT.
And OTFFT is a mixed-radix FFT.


[Caution]

OTFFT is optimized for the environment used to compile. In other words,
if we compile it on AVX-supported environment, it becomes the binary that uses
AVX. If someone executes this binary on AVX-unsupported environment, of course,
it crashes by causing the exception. For this reason, OTFFT is not suitable for
the use for distributing the compiled binary to all people. It is being assumed
that people use OTFFT to compile a program of numerical calculation from the
source code.


[Necessary files]

If you expand the "otfft-11.5e.tar.gz", you get several files.
Something necessary to you is "otfft" folder in these files.


[Setup]

In order to tune up this library, you need to execute the following commands in
"otfft" folder.

    make ffttune
    ./ffttune
    make otfft.o

For the compilation, it takes a very long time. Please be patient.
Some of linux will be very slow when it is run on 8 threads.
In that case, please set the "OMP_NUM_THREADS" environment variable to 7.


[How to use]

After you add the "otfft" folder in the compiler include path, please execute
the following commands. "hello.cpp" is your program.

    clang++-mp-4.0 -c -O3 hello.cpp
    clang++-mp-4.0 hello.o otfft/otfft.o -L/opt/local/lib/libomp -lomp -o hello

Here, I'm assuming that there is a "otfft" folder in the current folder.
"clang++-mp-4.0" is MacPotrs Clang compiler.


[Single Thread Mode]

To use the single thread mode, you need to compile the OTFFT with the condition
that DO_SINGLE_THREAD macro is defined at "otfft_misc.h". And, please start from
"make ffttune".

-------------------------------------------------------------------------------

Following such transforms are provided.

[Complex-to-Complex FFT]

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        complex_t* x = (complex_t*) simd_malloc(N*sizeof(complex_t));
        // do something
        OTFFT::FFT fft(N); // creation of FFT object. N is sequence length.
        fft.fwd(x);        // execution of transformation. x is input and output
        // do something
        simd_free(x);
    }

complex_t is defined as follows.

    struct complex_t
    {
        double Re, Im;

        complex_t() : Re(0), Im(0) {}
        complex_t(const double& x) : Re(x), Im(0) {}
        complex_t(const double& x, const double& y) : Re(x), Im(y) {}
        complex_t(const std::complex<double>& z) : Re(z.real()), Im(z.imag()) {}
        operator std::complex<double>() const { return std::complex<double>(Re, Im); }

        // ...
    };

There are member functions, such as the following.

    fwd(x)  -- DFT(with 1/N normalization)  x:input/output
    fwd0(x) -- DFT(non normalization)       x:input/output
    fwdu(x) -- DFT(unitary transformation)  x:input/output
    fwdn(x) -- DFT(with 1/N normalization)  x:input/output

    inv(x)  -- IDFT(non normalization)      x:input/output
    inv0(x) -- IDFT(non normalization)      x:input/output
    invu(x) -- IDFT(unitary transformation) x:input/output
    invn(x) -- IDFT(with 1/N normalization) x:input/output

To change the FFT size, do the following.

    fft.setup(2 * N);

To use in a multi-threaded environment, we do as follows.

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        complex_t* x = (complex_t*) simd_malloc(N*sizeof(complex_t));
        complex_t* y = (complex_t*) simd_malloc(N*sizeof(complex_t));
        // do someting
        OTFFT::FFT0 fft(N);
        fft.fwd(x, y); // x is input/output. y is work area
        // do something
        fft.inv(x, y); // x is input/output. y is work area
        // do someting
        simd_free(y);
        simd_free(x);
    }

Please note that "OTFFT::FFT" was changed to "OTFFT::FFT0".


[Real-to-Complex FFT]

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        double*    x = (double*)    simd_malloc(N*sizeof(double));
        complex_t* y = (complex_t*) simd_malloc(N*sizeof(complex_t));
        // do something
        OTFFT::RFFT rfft(N);
        rfft.fwd(x, y); // x is input. y is outout
        // do something
        simd_free(y);
        simd_free(x);
    }

N must be a multiple of 4. There are member functions, such as the following.

    fwd(x, y)  -- DFT(with 1/N normalization)  x:input, y:output
    fwd0(x, y) -- DFT(non normalization)       x:input, y:output
    fwdu(x, y) -- DFT(unitary transformation)  x:input, y:output
    fwdn(x, y) -- DFT(with 1/N normalization)  x:input, y:output

    inv(y, x)  -- IDFT(non normalization)      y:input, x:output
    inv0(y, x) -- IDFT(non normalization)      y:input, x:output
    invu(y, x) -- IDFT(unitary transformation) y:input, x:output
    invn(y, x) -- IDFT(with 1/N normalization) y:input, x:output

inv,inv0,invu,invn will destroy the input.


[Discrete Cosine Transformation(DCT-II)]

This transformation, orthogonalization is not executed.

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        double* x = (double*) simd_malloc(N*sizeof(double));
        // do something
        OTFFT::DCT dct(N);
        dct.fwd(x); // execution of DCT. x is input and output
        // do something
        simd_free(x);
    }

N must be an even number. There are member functions, such as the following.

    fwd(x)  -- DCT(with 1/N normalization) x:input/output
    fwd0(x) -- DCT(non normalization)      x:input/output
    fwdn(x) -- DCT(with 1/N normalization) x:input/output

    inv(x)  -- IDCT(non normalization)     x:input/output
    inv0(x) -- IDCT(non normalization)     x:input/output
    invn(x) -- IDCT(with 1/N nomalization) x:input/output

To use in a multi-threaded environment, we do as follows.

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        double*    x = (double*)    simd_malloc(N*sizeof(double));
        double*    y = (double*)    simd_malloc(N*sizeof(double));
        complex_t* z = (complex_t*) simd_malloc(N*sizeof(complex_t));
        // do something
        OTFFT::DCT0 dct(N);
        dct.fwd(x, y, z); // x is input/output. y,z are work area
        // do something
        dct.inv(x, y, z); // x is input/output. y,z are work area
        // do somthing
        simd_free(z);
        simd_free(y);
        simd_free(x);
    }

Please note that "OTFFT::DCT" was changed to "OTFFT::DCT0".


[Bluestein's FFT]

Bluestein's FFT is the FFT of any sequence length. Even if a sequence length is
a big prime number, the order of complexity is O(N log N).

    #include "otfft/otfft.h"
    using OTFFT::complex_t;
    using OTFFT::simd_malloc;
    using OTFFT::simd_free;

    void f(int N)
    {
        complex_t* x = (complex_t*) simd_malloc(N*sizeof(complex_t));
        // do something
        OTFFT::Bluestein bst(N);
        bst.fwd(x); // execution of Bluestein's FFT. x is input and output
        // do something
        simd_free(x);
    }

There are member functions, such as the following.

    fwd(x)  -- DFT(with 1/N normalization)  x:input/output
    fwd0(x) -- DFT(non normalization)       x:input/output
    fwdu(x) -- DFT(unitary transformation)  x:input/output
    fwdn(x) -- DFT(with 1/N normalization)  x:input/output

    inv(x)  -- IDFT(non normalization)      x:input/output
    inv0(x) -- IDFT(non normalization)      x:input/output
    invu(x) -- IDFT(unitary transformation) x:input/output
    invn(x) -- IDFT(with 1/N normalization) x:input/output

To use in a multi-threaded environment, we need to create objects of the same
number as the number of threads.

-------------------------------------------------------------------------------

[Benchmark]

Benchmark programs are bundled with OTFFT. These benchmarks require FFTW and
OOURA's FFT Package. Please refer to the following URL for OOURA's FFT Package.

    http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html

Please install the FFTW and place OOURA's "fftsg.c" file into "otfft-11.5e"
folder. Then, execute the following command.

    make fftbench1

And, execute the following command to execute the benchmark.

    ./fftbench1
