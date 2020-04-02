

## Build Instructions
```
#Prerequisites:
$ sudo apt install ftjam libglm-dev libsdl2-dev libpng-dev

#Checkout code + submodules:
$ git clone git@github.com:ixchow/freiscale
$ cd freiscale
frescale$ git submodule update --init

#Tune otfft (optional):
freiscale$ cd otfft-11.5e/otfft
freiscale/otfft-11.5e/otfft$ make ffttune
freiscale/otfft-11.5e/otfft$ ./ffttune
freiscale/otfft-11.5e/otfft$ cd ../..

#Build otfft (not optional):
freiscale$ cd otfft-11.5e/otfft
freiscale/otfft-11.5e/otfft$ make otfft.o
freiscale/otfft-11.5e/otfft$ cd ../..

#Build:
freiscale$ jam
```

You can also build against the `ixchow/kit-libs-linux`, `ixchow/kit-libs-osx`, or `ixchow/kit-libs-win` packages instead of installing glm, sdl2, and libpng system-wide. This is especially useful on windows. Just check out the corresponding repository as a subdirectory of freiscale.
