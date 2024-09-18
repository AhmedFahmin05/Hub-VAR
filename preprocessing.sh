#!/bin/bash
f="$(basename -- $1)";
distance=$1"/USA-road-d."$f".gr";
converted=$1"/"$f".txt";
./bin/convert_graph $distance $converted
./bin/convert $distance $converted


output=$1"/results/"
label=$1"/"$f;

./bin/construct $converted $label $output
