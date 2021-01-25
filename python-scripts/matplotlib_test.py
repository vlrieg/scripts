#!/usr/bin/env python3


import matplotlib.pyplot as plt

x = [i for i in range(101)]
y = [i**2 for i in x]

a = [i for i in range(51)]
b = [i*2 for i in a]

my_vals = [
    (x, y),
    (a, b)
]

index = 0
for vals in my_vals:
    plt.plot(vals[0], vals[1], label=f"value {index}")
    index += 1

plt.legend()
plt.savefig("matplotlib_test_output.pdf")

