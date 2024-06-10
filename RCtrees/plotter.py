import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run RC program and plot results.')
parser.add_argument('--randomized', type=str, choices=['true', 'false'], default='false', help='Set randomized flag for RC.')
parser.add_argument('--do-height', type=str, choices=['true', 'false'], default='true', help='Set do-height flag for RC.')
parser.add_argument('-n', type=str, required=True, help='Name to append to titles and file names')
args = parser.parse_args()

# Constants
graph_sizes = [2500000000, 1000000000, 500000000, 100000000, 50000000, 10000000, 5000000, 1000000, 100000, 10000, 1000, 100]


num_threads_list = [144, 100, 72, 60, 48, 36, 24, 16, 8, 4, 2, 1]
name_suffix = args.n
time_output_filename = f'creation_time_vs_graph_size_randomized_{args.randomized}_do_height_{args.do_height}_{name_suffix}.png'
speedup_output_filename = f'speedup_vs_graph_size_randomized_{args.randomized}_do_height_{args.do_height}_{name_suffix}.png'

# Prepare to collect results
results = {threads: [] for threads in num_threads_list}

# Function to run the C++ program with a given number of threads and graph size
def run_rc_with_threads(graph_size, num_threads):
    env = os.environ.copy()
    env['PARLAY_NUM_THREADS'] = str(num_threads)
    env['LD_PRELOAD'] = '/usr/local/lib/libjemalloc.so'
    try:
        cmd = ['./RC', '--print-creation', '--randomized', args.randomized, '--do-height', args.do_height, '-n', str(graph_size)]
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        output = result.stdout.strip()
        if output:
            size, time = output.split(',')
            return int(size), float(time)
    except Exception as e:
        print(f"Failed to run ./RC for graph size {graph_size} with {num_threads} threads: {e}")
    return None

# Collect results for each configuration
for num_threads in num_threads_list:
    for graph_size in graph_sizes:
        print(f'{graph_size=}')
        result = run_rc_with_threads(graph_size, num_threads)
        if result:
            results[num_threads].append(result)
            print(f"Graph size: {result[0]}, Time: {result[1]} seconds with {num_threads} threads")

# Prepare data for plotting
plot_data = {threads: {'sizes': [], 'times': []} for threads in num_threads_list}
for num_threads, data in results.items():
    for size, time in data:
        plot_data[num_threads]['sizes'].append(size)
        plot_data[num_threads]['times'].append(time)

# Calculate speedup relative to single-threaded execution
speedup_data = {threads: [] for threads in num_threads_list if threads != 1}
single_thread_times = {size: time for size, time in results[1]}

for num_threads, data in results.items():
    if num_threads == 1:
        continue
    for size, time in data:
        if size in single_thread_times:
            speedup_data[num_threads].append(single_thread_times[size] / time)

# Plot creation times vs graph sizes
plt.figure(figsize=(12, 8))
for num_threads in num_threads_list:
    plt.plot(plot_data[num_threads]['sizes'], plot_data[num_threads]['times'], marker='o', label=f'{num_threads} threads')
plt.xlabel('Graph Size (n)')
plt.ylabel('Creation Time (seconds)')
plt.title(f'Graph Creation Time vs Graph Size (randomized={args.randomized}, do-height={args.do_height}) {name_suffix}')
plt.legend()
plt.savefig(time_output_filename)
print(f'Creation time plot saved as {time_output_filename}')

# Plot speedup vs graph sizes
plt.figure(figsize=(12, 8))

# Generate constant spacing for the x-axis
x_positions = np.arange(len(num_threads_list) - 1)  # Exclude single thread from speedup plot

# Plot speedup with constant spacing
for i, num_threads in enumerate(num_threads_list):
    if num_threads != 1:
        plt.plot(plot_data[num_threads]['sizes'], speedup_data[num_threads], marker='x', linestyle='--', label=f'Speedup {num_threads} threads')

plt.xlabel('Graph size')
plt.xscale('log')
plt.ylabel(f'Speedup (relative to 1 thread) (randomized={args.randomized}, do-height={args.do_height}) {name_suffix}')
plt.title(f'Speedup vs Graph Size (randomized={args.randomized}, do-height={args.do_height}) {name_suffix}')
plt.legend()
plt.savefig(speedup_output_filename)
print(f'Speedup plot saved as {speedup_output_filename}')
