# Need to have numpy, pandas and matplotlib installed for this script to work.
# aka: `pip install pandas numpy matplotlib`

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Not currently the neatest of files..

# methods = [ "hpce.fast_fourier_transform", "hpce.direct_fourier_transform", "hpce.txl11.fast_fourier_transform_taskgroup", 
# "hpce.txl11.fast_fourier_transform_parfor", "hpce.txl11.fast_fourier_transform_combined", "hpce.txl11.fast_fourier_transform_opt" ]

# legend = []

# for m in methods:
#     name = "final/{0}.csv".format(m)

#     df=pd.read_csv(name, sep=',',header=None)

#     # print df.values
#     times = df.values[:,5][1:].astype(np.float_)
#     ns = df.values[:,3][1:].astype(np.float_)

#     plt.plot(ns, times)
#     legend.append(m)

# plt.yscale('log')
# plt.xscale('log')

# plt.title('Different FFT Methods')
# plt.xlabel('n (data length)')
# plt.ylabel('Execution time (seconds)')

# plt.legend(legend, loc='upper right')
# plt.show()


processors = range(1,9)
# Ks = [np.power(2, x) for x in range(0, 13)]
# Ks = [2, 64, 128, 256, 1024, 4096]
Ks = [512]

legend = []

for k in Ks:
    for procs in processors:
        name = "fft_taskgroup_K_{0}_P_{1}.csv".format(k,procs)

        df=pd.read_csv(name, sep=',',header=None)

        # print df.values
        times = df.values[:,5][1:].astype(np.float_)
        ns = df.values[:,3][1:].astype(np.float_)

        plt.plot(ns, times)
        legend.append("P={0}".format(procs))

plt.yscale('log')
plt.xscale('log')

plt.title('FFT Taskgroup - Number of Processors only (P). All K=512')
plt.xlabel('n (data length)')
plt.ylabel('Execution time (seconds)')

plt.legend(legend, loc='upper right')
plt.show()



# processors = ['1', '3', '8']
# # Ks = [np.power(2, x) for x in range(0, 13)]
# Ks = [2, 64, 128, 256, 1024, 4096]

# legend = []

# for k in Ks:
#     for procs in processors:
#         name = "fft_taskgroup_K_{0}_P_{1}.csv".format(k, procs)

#         df=pd.read_csv(name, sep=',',header=None)

#         # print df.values
#         times = df.values[:,5][1:].astype(np.float_)
#         ns = df.values[:,3][1:].astype(np.float_)

#         plt.plot(ns, times)
#         legend.append("K={0}, P={1}".format(k, procs))

# plt.yscale('log')
# plt.xscale('log')

# plt.title('FFT Taskgroup - Varying K and Number of Processors (P)')
# plt.xlabel('n (data length)')
# plt.ylabel('Execution time (seconds)')

# plt.legend(legend, loc='upper right')
# plt.show()


# processors = ['1', '3', '8']
# # Ks = [np.power(2, x) for x in range(0, 13)]
# Ks = [2, 64, 128, 256, 1024, 4096]

# legend = []

# for k in Ks:
#     for procs in processors:
#         name = "fft_parfor_K_{0}_P_{1}.csv".format(k, procs)

#         df=pd.read_csv(name, sep=',',header=None)

#         # print df.values
#         times = df.values[:,5][1:].astype(np.float_)
#         ns = df.values[:,3][1:].astype(np.float_)

#         plt.plot(ns, times)
#         legend.append("K={0}, P={1}".format(k, procs))

# plt.yscale('log')
# plt.xscale('log')

# plt.title('FFT Parfor (Chunking) - Varying K and Number of Processors (P)')
# plt.xlabel('n (data length)')
# plt.ylabel('Execution time (seconds)')

# plt.legend(legend, loc='upper right')
# plt.show()