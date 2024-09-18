#!/bin/bash
f="$(basename -- $1)";
distance=$1"/USA-road-t."$f".gr";
converted=$1"/"$f"_time";
./bin/convert_graph $distance $converted
./bin/converter $distance $converted


output=$1"/results/"
label=$1"/"$f;

./bin/construct $converted $label $output
