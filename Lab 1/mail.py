#%%
import matplotlib.pyplot as plt

#%%
import numpy as np
#%%
with open('f14.txt') as file:
    data = np.array([float(val) for val in file.read().split()])

T = 5
dt = 0.01
time = np.arange(0, T + dt, dt)
plt.grid(True)
plt.plot(time, data)
#%%
