#!/bin/sh -

for f in *.h *.cpp LICENSE.txt
do
    mv $f $f.bak && nkf --windows $f.bak > $f && rm $f.bak
done
