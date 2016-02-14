#estimating the time complexity of the thermal modeling using the time module

import matplotlib.pyplot as plt
import time
import numpy as np
from thermal_model import run_simulation 

num_times = 30
days = 0.03 #days

error_space = np.linspace(10**(-8), 10**(-4), num=num_times) 
results = []
for rtol in error_space:
    then = time.time()
    run_simulation(days, rtol)
    now = time.time()
    results.append(now - then)
    print("simulation completed for days: " + str(days))
    print("Time: " + str(now - then))
    
plt.plot(error_space, results)
plt.xlabel("time in days")
plt.ylabel("time elapsed for simulation")
plt.show()