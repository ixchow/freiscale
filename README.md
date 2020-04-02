# Freiscale Semicomposer

A (very rough prototype of a) composition program with extensive visualization and no scales.

## Suggested Sample Library

If you launch `freiscale` with the command-line parameter `lib:path/` it will look in `path/` (relative to the current working directory) for samples. If you don't specify the command-line parameter, it will look in the directory `sounds/` (relative to the executable) for a sample library.

A sample library is just a bunch of `wav` files in folders. You might want to give them descriptive names.

The OLPC project has links to samples here:
http://wiki.laptop.org/go/Free_sound_samples

Particularly, I have used the "Open Path Music Collection" packages in much of my testing.

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

NOTE: If you are getting weird segfaults in the FFT code, there might be a compile flags mis-match between the flags used for otfft.o and the flags used by the code that interfaces with otfft.o .
