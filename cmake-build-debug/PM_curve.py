import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
import binascii

hist_points = [5, 50, 500]
f = open('snapshots/exper_512_configs.txt', 'r').readlines()
# ar_sum = sum([np.exp(-1*int(x)) for x in f])
runs = [int(x, 2) for x in f]
#print(f)
x1, x2, x3 = runs[4::500], runs[49::500], runs[499::500]
# hist_steps = [[x[hist_points[0]-1] for x in runs], [x[hist_points[1]-1] for x in runs], [x[hist_points[2]-1] for x in runs]]
for x in [x1, x2, x3]:
    plt.hist(x)
    plt.title('Experimental P(M)')
    plt.xlabel('Configuration')
    plt.ylabel('Probability')
    plt.show()

# plt.errorbar(x=range(0, len(runs[20:]), 5),y=runs[20::5], yerr=[.001 for x in range(0, len(runs[20::5]))], fmt=',')

# for step in hist_steps:
#     plt.hist(step)
#     plt.show()

#
# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
