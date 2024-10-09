import numpy as np
import matplotlib.pyplot as plt

def theoretical_arithmetic_intensity(nnz, nnode, nblock):
    # byte = (4 + 8 * (9 + 3*nblock)) * nnz + 8 * 3*nblock * nnode
    byte = (4 + 8 * (9 + 3)) * nnz + 8 * 3 * nnode
    flop = 18 * nnz * nblock
    return flop / byte

nnode, nnz, = 2146689, 5706662

# Define the parameters for the roofline model
peak_flops = 7.362e12  # Peak FLOPs/s
peak_bandwidth = 1935e9  # Peak memory bandwidth in bytes/s

# Define the operational intensity range
operational_intensity = np.logspace(-2, 2, 100)

# Calculate the performance for each operational intensity
performance = np.minimum(peak_flops, operational_intensity * peak_bandwidth)

# Plot the roofline model
plt.figure(figsize=(10, 6))
plt.loglog(operational_intensity, performance, label='Roofline', linewidth=2)

# Add labels and title
plt.xlabel('Operational Intensity (FLOPs/Byte)')
plt.ylabel('Performance (FLOPs/s)')
plt.title('GPU Roofline Model')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# from sayson's thesis
plt.scatter([0.19, 0.38, 0.72], [0.0224*9.7e12, 0.0466*9.7e12, 0.0855*9.7e12], color='red', label='SGEMM')

# result
plt.scatter([0.19, 0.40, 0.77, 1.35, 2.15], [0.225e12, 0.447e12, 0.817e12, 1.264e12, 1.665e12], color='blue', label='SGEMM')

# Show the plot
plt.savefig("Roofline.pdf")