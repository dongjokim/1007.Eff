#!/bin/bash

GREPARG=${1:-".*"}

find data | grep LHC  | grep "$GREPARG" |  while read DIR;
do
    #TAG=`echo $DIR | perl -nle'/__([a-zA-Z0-9]+)/ and print $1 or /-(r\w+)/ and print $1'`
    TAG=`echo $DIR | perl -nle'/_(.+?)\./ and print $1 or /-(r\w+)/ and print $1'`
    echo root -b -q -l calEffTrigg.C+g\'\(\"$DIR\",\"$TAG\"\)\'
done

