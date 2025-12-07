import csv
from matplotlib import pyplot as plt
import numpy as np
import math

x = []
y = []

with open('build/numerical.csv', mode ='r')as file:
    csv_file = csv.reader(file)

    for line in csv_file:
        x.append(float(line[0]))
        y.append(float(line[1]))

y0 = 2
analytic_x = np.arange(0,5,0.005)
analytic_y = []
for xs in analytic_x:
    analytic_y.append(math.exp(0.5*(xs-0.5*math.sin(2*xs))+math.log(y0)))

plt.plot(x,y, "x-")
plt.plot(analytic_x, analytic_y)
plt.show()