#plot random vs deterministic

import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np

# Constants
graph_size = 100000000
num_threads_list = [1, 2] + [i for i in range(4, 145, 4)]
print(f'testing for {num_threads_list=}')
randomized_values = ['true', 'false']
do_height = 'false'
time_output_filename = 'creation_time_vs_num_threads max_deg=3.png'
speedup_output_filename = 'speedup_vs_num_threads max_deg=3.png'

# Function to run the C++ program with a given number of threads and randomized setting
def run_rc_with_threads(graph_size, num_threads, randomized):
    env = os.environ.copy()
    env['PARLAY_NUM_THREADS'] = str(num_threads)
    env['LD_PRELOAD'] = '/usr/local/lib/libjemalloc.so'
    try:
        cmd = ['./RC', '--print-creation', '--randomized', randomized, '--do-height', do_height, '-n', str(graph_size)]
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        output = result.stdout.strip()
        if output:
            size, time = output.split(',')
            return float(time)
    except Exception as e:
        print(f"Failed to run ./RC for {num_threads} threads (randomized={randomized}): {e}")
    return None


# Collect results for each configuration
results = {randomized: {threads: [] for threads in num_threads_list} for randomized in randomized_values}
for randomized in randomized_values:
    for num_threads in num_threads_list:
        print(f'{num_threads=}')
        result = run_rc_with_threads(graph_size, num_threads, randomized)
        if result:
            results[randomized][num_threads] = result
            print(f"Threads: {num_threads}, Randomized: {randomized}, Time: {result} seconds")

# Plot creation times vs number of threads
plt.figure(figsize=(12, 8))
for randomized in randomized_values:
    plt.plot(num_threads_list, [results[randomized][threads] for threads in num_threads_list], marker='o', label=f'Randomized={randomized}')
plt.xlabel('Number of Threads')
plt.ylabel('Creation Time (seconds)')
plt.title(f'Creation Time vs Number of Threads (Graph Size: {graph_size}, do-height={do_height}) max_deg=3')
plt.legend()
plt.savefig(time_output_filename)
print(f'Creation time plot saved as {time_output_filename}')

# Plot speedups vs number of threads
plt.figure(figsize=(12, 8))
for randomized in randomized_values:
    speedup_data = [results[randomized][1] / results[randomized][threads] for threads in num_threads_list]
    plt.plot(num_threads_list, speedup_data, marker='x', linestyle='--', label=f'Speedup (Randomized={randomized})')
plt.xlabel('Number of Threads')
plt.ylabel('Speedup (relative to 1 thread)')
plt.title(f'Speedup vs Number of Threads (Graph Size: {graph_size}, do-height={do_height}) max_deg=3')
plt.legend()
plt.savefig(speedup_output_filename)
print(f'Speedup plot saved as {speedup_output_filename}')
