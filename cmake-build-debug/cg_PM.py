import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
import binascii

f = open('snapshots/PM_cg1.txt', 'r').readlines()
# ar_sum = sum([np.exp(-1*int(x)) for x in f])
M = [int(x) for x in f[::2]]
configs = [int(x, 2) for x in f[1::2]]
# for i in range(0, len(M)):
#     PMs = [(x[0], x[1]) for x in zip(M, configs)]
#     l = [[x[0] for x in PMs if x[0] == 0], [x[0] for x in PMs if x[0] == 2], [x[0] for x in PMs if x[0] == 4]]
print(M[:5], configs[:5])
plt.hist(configs)
plt.title('Experimental P(M) N/3')
plt.xlabel('Configuration')
plt.ylabel('Probability')
plt.show()
