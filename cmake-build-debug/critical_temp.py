import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
import binascii

hist_points = [5, 50, 500]
steps = 5
f = open('snapshots/critical_temp_3.txt', 'r').readlines()
f += open('snapshots/critical_temp_5.txt', 'r').readlines()
f += open('snapshots/critical_temp_9.txt', 'r').readlines()
f += open('snapshots/critical_temp_15.txt', 'r').readlines()


# ar_sum = sum([np.exp(-1*int(x)) for x in f])
m = [float(x) for x in f[::2]]
betas = [float(x) for x in f[1::2]]

for i in range(0, len(m), steps):
    plt.plot(betas[i:i+steps], m[i:i+steps])

plt.title('< m^2 > vs β')
plt.xlabel('β')
plt.ylabel('< m^2 >')
plt.show()

# plt.errorbar(x=range(0, len(runs[20:]), 5),y=runs[20::5], yerr=[.001 for x in range(0, len(runs[20::5]))], fmt=',')

# for step in hist_steps:
#     plt.hist(step)
#     plt.show()

#
# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
