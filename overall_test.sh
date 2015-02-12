#!/bin/bash

TASKS="hpce.fast_fourier_transform hpce.direct_fourier_transform"
TASKS="${TASKS} hpce.txl11.fast_fourier_transform_taskgroup" 
TASKS="${TASKS} hpce.txl11.fast_fourier_transform_parfor hpce.txl11.fast_fourier_transform_combined hpce.txl11.fast_fourier_transform_opt"


for t in $TASKS; do
    echo "Running test for ${t}";
    ./bin/time_fourier_transform  $t > test/final/$t.csv
done