#!/bin/sh -

for f in *.h *.cpp *.f03 LICENSE.txt
do
    mv $f $f.bak && nkf --windows $f.bak > $f && rm $f.bak
done
