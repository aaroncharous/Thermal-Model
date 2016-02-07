#estimating the time complexity of the thermal modeling using the time module

import matplotlib.pyplot as plt
import time
import numpy as np
from thermal_model import run_simulation 

num_times = 10
max_time = 0.25 #days

time_space = np.linspace(0.025, max_time, num=num_times) 
results = []
for days in time_space:
    then = time.time()
    run_simulation(days, graphing=False)
    now = time.time()
    results.append(now - then)
    print("simulation completed")
    
plt.plot(time_space, results)
xlabel("time in days")
ylabel("time elapsed for simulation")
#plt.show()    





