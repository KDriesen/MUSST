#-*- coding: utf-8 -*-

import sys
import math
import numpy as np
from  matplotlib import pyplot as plt

#read data from a file
#name of the variables
#unpack is true = variables in columns
x, h, p, pt = np.loadtxt(sys.argv[1], unpack = True)
#plot of the data
plt.plot(x, p, marker = 'o', linestyle='None')
plt.plot(x, pt)
#plt.xlabel('$x/L$')
#plt.ylabel('$ph^2/6 \mu V L$')
plt.xlabel(sys.argv[2])
plt.ylabel(sys.argv[3])
plt.show()
