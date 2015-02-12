#!/bin/bash

TASKS="hpce.fast_fourier_transform hpce.direct_fourier_transform hpce.txl11.direct_fourier_transform_parfor"
TASKS="${TASKS} hpce.txl11.direct_fourier_transform_parfor hpce.txl11.fast_fourier_transform_taskgroup hpce.txl11.direct_fourier_transform_chunked" 
TASKS="${TASKS} hpce.txl11.fast_fourier_transform_parfor hpce.txl11.fast_fourier_transform_combined hpce.txl11.fast_fourier_transform_opt"

for t in $TASKS; do
    echo "Running test for ${t}";
    ./bin/test_fourier_transform $t
done
