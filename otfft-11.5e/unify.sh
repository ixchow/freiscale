#!/bin/sh -

for f in otfft/otfft_gen_*.h
do
    mv $f $f.tmp && sed "s/FFT[1-7]/FFT$1/" $f.tmp > $f && rm $f.tmp
done
