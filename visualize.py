import csv
from matplotlib import pyplot as plt
import numpy as np
import math

x = []
y = []
y2 = []

with open('build/numerical.csv', mode ='r')as file:
    csv_file = csv.reader(file)

    for line in csv_file:
        x.append(float(line[0]))
        y.append(float(line[1]))
        y2.append(math.pow(float(line[1]),2))

y0 = 2
analytic_x = np.arange(0,5,0.005)
analytic_y = []
for xs in analytic_x:
    analytic_y.append(math.exp(0.5*(xs-0.5*math.sin(2*xs))+math.log(y0)))

plt.plot(x,y, "-")
plt.plot(x,y2)
plt.xlabel("Distance (m)")
plt.ylabel("Value")
plt.legend(["Wavefunction", "Probability"])
#plt.plot(analytic_x, analytic_y)
plt.show()

# Params
params = []
errors = []

with open('build/params.csv', mode ='r')as file:
    csv_file = csv.reader(file)

    for line in csv_file:
        params.append(float(line[0]))
        errors.append(float(line[1]))

plt.plot(params, errors, "x-")
plt.xlabel("Param")
plt.ylabel("Error")
plt.legend(["Error"])
#plt.plot(analytic_x, analytic_y)
plt.show()
