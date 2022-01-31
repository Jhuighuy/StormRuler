#!/bin/bash

awk 'BEGIN { i = 1 } $1 == "V_hat" { print i, $3, $6; i += 1 }' slices.txt > slices_plot.txt

#first is computed, second -- reference..
#plot 'slices_plot.txt' u 1:2, 'slices_plot.txt' u 1:3