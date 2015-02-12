#!/bin/bash

# This assumes an 8 core machine..
CORES="1 2 3 4 5 6 7 8";

K="1 2 4 8 16 32 64 128 256 512 1024 2048 4096";

# direct_fourier_transform_chunked
for procs in $CORES; do
    for val in $K; do
        # echo $K;
        export HPCE_DIRECT_OUTER_K=$val;
        fname="ft_chunked_K_${val}_P_${procs}"
        ./bin/time_fourier_transform hpce.txl11.direct_fourier_transform_chunked $procs 60 > test/$fname.csv 2> test/error/$fname.log
    done
done


# fast_fourier_transform_taskgroup
for procs in $CORES; do
    for val in $K; do
        fname="fft_taskgroup_K_${val}_P_${procs}"
        export HPCE_FFT_RECURSION_K=$val;
        ./bin/time_fourier_transform hpce.txl11.fast_fourier_transform_taskgroup $procs 60 > test/$fname.csv 2> test/error/$fname.log
    done
done


# fast_fourier_transform_parfor
for procs in $CORES; do
    for val in $K; do
        fname="fft_parfor_K_${val}_P_${procs}"
        export HPCE_FFT_LOOP_K=$val;
        ./bin/time_fourier_transform hpce.txl11.fast_fourier_transform_parfor $procs 60 > test/$fname.csv 2> test/error/$fname.log
    done
done


# fast_fourier_transform_combined
for procs in $CORES; do
    for val in $K; do
        for innerval in $K; do
            fname="fft_combined_REC_K_${val}_LOOP_K_${innerval}_P_${procs}"
            export HPCE_FFT_RECURSION_K=$val;
            export HPCE_FFT_LOOP_K=$innerval;
            ./bin/time_fourier_transform hpce.txl11.fast_fourier_transform_parfor $procs 60 > test/fname.csv 2> test/error/fname.log
        done
    done
done