# coding: utf-8

# In[1]:


import numpy as np
import csv
import matplotlib.pyplot as plt


# In[2]:


def gaussian(pulse_width, center, time):
    sigma = pulse_width / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-(time - center) ** 2 / (2 * sigma ** 2))


# ### Sampling parameters. Do not change unless planning to change sampling rate

# In[3]:


# setting AWG params
sampling_rate = 92e9
sample_time = 1 / sampling_rate

# ### Change the pulse width , seperation and overall pattern repetition time

# In[4]:


# setting pulse times
pulse_width = 200e-12
pulse_seperation = 2e-9
pulse_shape = "gaussian"
pulse_repetition = 500e-9

# In[13]:


total_data_points = int((pulse_width + pulse_seperation + pulse_repetition) // sample_time)

if total_data_points > 512000:
    print("Warning. The waveform is bigger than AWG's memory")

# In[15]:


file_name = "AWG" + "_timebin_" + str(pulse_width // 1e-12) + "ps" + "_sep_" + str(
    pulse_seperation // 1e-9) + "ns" + "_" + pulse_shape

# In[22]:


pulse_width_in_steps = pulse_width // sample_time
pulse_seperation_in_steps = pulse_seperation // sample_time
center = 100
x = []
y = []
m = []
p_first = []
p_second = []
for i in range(total_data_points):
    x.append(i)
    amp = gaussian(pulse_width_in_steps, center, i) + gaussian(pulse_width_in_steps,
                                                               (center + pulse_seperation_in_steps), i)
    if amp > 5e-5:
        y.append(amp)
        if i < (center + (pulse_width_in_steps * 2)):
            p_first.append(1)
            p_second.append(0)

        else:
            p_first.append(0)
            p_second.append(1)

    else:
        y.append(0)
        p_first.append(0)
        p_second.append(0)
    m.append(0)

# In[23]:


fig, ax = plt.subplots(figsize=(20, 5))
ax.plot(x, y, x, p_first, 'r--', x, p_second, 'g--')
plt.xlim(0, 500)
fig.savefig(file_name + '.png')
plt.show()

# In[121]:


with open(file_name + '.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(y, m))

# ##Outputs phase modulating square pulses. *phase_1* in file name notes square pulse on first peak and *phase_2* in file name notes square pulse on second peak

# In[27]:


with open(file_name + '_phase_1' + '.csv', 'w') as pp:
    writer = csv.writer(pp)
    writer.writerows(zip(p_first, m))
with open(file_name + '_phase_2' + '.csv', 'w') as pm:
    writer = csv.writer(pm)
    writer.writerows(zip(p_second, m))

# In[ ]:

