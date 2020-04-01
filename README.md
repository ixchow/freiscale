* >> suo apt install jam
* >> sudo apt install libglm-dev
* >> sudo apt install libsdl2-dev
* >> sudo apt install libpng-dev
* >> sudo apt install libfftw3-dev
* clone the freiscale repo, enter the directory
* >> git submodule update --init
* >> cd otfft-11.5e
* >> wget http://www.kurims.kyoto-u.ac.jp/~ooura/fft.tgz
* >> tar -xf fft.tgz
* >> mv fft/fftsg.c .
* >> cd ../
* >> jam