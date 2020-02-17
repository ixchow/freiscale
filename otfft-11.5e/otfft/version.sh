#!/bin/sh -

for f in otfft*.h otfft.cpp otfft_c.cpp otfft.f03 ffttune.cpp rewrite.cpp
do
    mv $f $f.tmp && \
        sed "2s/Version [1-9][0-9]*[.0-9]*[a-z]*/Version $1/" $f.tmp > $f && \
        rm $f.tmp
done
