
import numpy as np
import matplotlib.pyplot as plt

def phi(s1, s2):
    return np.degrees(np.arccos((s1**2/2 + s2**2/2 + 2*s1*s2) / (s1**2 + s2**2 + s1*s2)))

s1_values = np.linspace(-100, 100, 1000)  
s2_values = np.arange(-11, 11, 1)

plt.figure(figsize=(10, 6))

for s2 in s2_values:
    phi_values = phi(s1_values, s2)
    plt.plot(s1_values, phi_values, label=f's2 = {s2}')

plt.axhline(y=1.08, color='r', linestyle='--', label='Magic Angle: 1.08°')

plt.xlabel('s1')
plt.ylabel('phi (degrees)')
plt.title('Plot of phi as a function of s1 for different s2 values')
plt.legend()
plt.grid()
plt.show()

