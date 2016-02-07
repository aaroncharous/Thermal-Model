
# coding: utf-8

# In[12]:

#estimating the time complexity of the thermal modeling using the time module

import matplotlib.pyplot as plt
import time
import numpy as np
from thermal_model import run_simulation 

num_times = 10
max_time = 0.25 #days

time_space = np.linspace(0, max_times, num=num_times) 
results = []
print(time_space)
for days in time_space:
    then = time.time()
    run_simulation(days, graph=False)
    now = time.time()
    results.append(now - then)
    
plt.plot(time_space, results)
xlabel("time in days")
ylabel("time elapsed for simulation")
    



# In[10]:

import time


# In[13]:

import matplotlib.pyplot as plt


# 

# In[ ]:



