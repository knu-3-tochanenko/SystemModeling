# %%
# Import libraries
import numpy as numpy
import matplotlib.pyplot as plot
# %%
# Open file
with open('f14.txt') as file:
    inputData = numpy.array([float(val) for val in file.read().split()])
# %%
# Build det plot
T = 5
dt = 0.01
time = numpy.arange(0, T+dt, dt)
plot.grid(True)
plot.rc('figure', figsize=(30.0, 10.0))
plot.plot(time, inputData)
# %%
# Build time plot
n = time.shape[0]
frequency = numpy.zeros(n)
for pointID in range(n):
    sinFrequency = 0
    cosFrequency = 0
    for signal in range(n):
        sinFrequency += inputData[signal] * \
            numpy.sin(2.*numpy.pi*pointID*signal/float(n))
        cosFrequency += inputData[signal] * \
            numpy.cos(2.*numpy.pi*pointID*signal/float(n))
    sinFrequency /= float(n)
    cosFrequency /= float(n)
    frequency[pointID] = numpy.sqrt(sinFrequency**2+cosFrequency**2)
plot.grid(True)
plot.plot(time, frequency)

# %%
# Find biggest value
biggestValue = []
for i in range(3, n // 2):
    if numpy.max(frequency[i-3:i+3]) == frequency[i]:
        biggestValue.append(i)
        print(frequency[i])
mainFrequency = biggestValue[0]/T
