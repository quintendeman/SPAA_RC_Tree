import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV data into a pandas DataFrame
data = pd.read_csv('output.csv')

# Find the global maximum of the 'num_threads' column
max_threads = data['num_threads'].max()

# Filter the data to include only rows where 'num_threads' is the global maximum
filtered_data = data[data['num_threads'] == max_threads]

# Group the filtered data by 'ln', 'mean', 'dist', and 'randomize'
grouped_data = filtered_data.groupby(['ln', 'mean', 'dist', 'randomize'])

# Create a figure for the static time plot
plt.figure(figsize=(10, 6))

# Loop through each group and plot for static_gen_time and static_tern_time
for (ln, mean, dist, randomize), group in grouped_data:
    # Replace 'e' with 'Exponential' and 'u' with 'Uniform' in 'dist'
    if dist == "e":
        dist = "Exponential"
    elif dist == "u":
        dist = "Uniform"
    
    # Sort the group by 'graph_size'
    group = group.sort_values('graph_size')
    
    # Generate a unique color for each line
    color = np.random.rand(3,)  # Random RGB color
    
    # Plot static_gen_time (non-dotted line)
    plt.plot(group['graph_size'], group['static_gen_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color)
    
    # Plot static_tern_time (dotted line)
    plt.plot(group['graph_size'], group['static_tern_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color, linestyle='dotted')

# Customize the static plot
plt.title('Static Generation and Ternerization Time vs Graph Size (Max Threads)')
plt.xlabel('Graph Size')
plt.ylabel('Time (s)')
plt.grid(True)

# Adjust the legend to be outside, to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Configurations", fontsize=8)

# Save the static time figure
plt.savefig('graphs/static_time.png', bbox_inches='tight')

# Show the static plot
plt.show()

# Create a figure for the dynamic time plot
plt.figure(figsize=(10, 6))

# Loop through each group and plot for dynamic_gen_time and dynamic_tern_time
for (ln, mean, dist, randomize), group in grouped_data:
    # Replace 'e' with 'Exponential' and 'u' with 'Uniform' in 'dist'
    if dist == "e":
        dist = "Exponential"
    elif dist == "u":
        dist = "Uniform"
    
    # Sort the group by 'graph_size'
    group = group.sort_values('graph_size')
    
    # Calculate batch-insert-size as graph_size / mean
    batch_insert_size = group['graph_size'] / group['mean']
    
    # Generate a unique color for each line
    color = np.random.rand(3,)  # Random RGB color
    
    # Plot dynamic_gen_time (non-dotted line)
    plt.plot(batch_insert_size, group['dynamic_gen_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color)
    
    # Plot dynamic_tern_time (dotted line)
    plt.plot(batch_insert_size, group['dynamic_tern_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color, linestyle='dotted')

# Customize the dynamic plot
plt.title('Dynamic Insertiong and Ternerization Time vs Batch Insert Size (Max Threads)')
plt.xlabel('Batch Insert Size (graph_size / mean)')
plt.ylabel('Time (s)')
plt.grid(True)

# Adjust the legend to be outside, to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Configurations", fontsize=8)

# Save the dynamic time figure
plt.savefig('graphs/dynamic_time.png', bbox_inches='tight')

# Show the dynamic plot
plt.show()

# Group data by 'ln', 'mean', 'dist', and 'randomize'
grouped_data = data.groupby(['ln', 'mean', 'dist', 'randomize'])

# Create a figure for the static time plot (for static_gen_time vs num_threads)
plt.figure(figsize=(10, 6))

unique_graph_size = 0
unique_mean = 0

# Loop through each group
for (ln, mean, dist, randomize), group in grouped_data:
    # Check if the configuration has multiple num_threads
    if len(group['num_threads'].unique()) > 1:
        # Filter the group to use static_gen_time and plot against num_threads
        group = group.sort_values('num_threads')
        
        # Generate a unique color for each line
        color = np.random.rand(3,)  # Random RGB color

        unique_graph_size = np.median(group['graph_size'])
        
        # Plot static_gen_time against num_threads
        plt.plot(group['num_threads'], group['static_gen_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color)

# Customize the static plot
plt.title(f'Static Generation Time vs Number of Threads (Graph Size={unique_graph_size})')
plt.xlabel('Number of Threads')
plt.ylabel('Static Generation Time (s)')
plt.grid(True)

# Adjust the legend to be outside, to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Configurations", fontsize=8)

# Save the static time figure
plt.savefig('graphs/static_threads_time.png', bbox_inches='tight')

# Show the static plot
plt.show()

# Find the median of graph_size/mean
median_batch_insert_size = data['graph_size'] / data['mean']
median_batch_insert_size = median_batch_insert_size.median()

# Filter the data for the median batch insert size
filtered_data_dynamic = data[median_batch_insert_size == (data['graph_size'] / data['mean'])]

# Group data by 'ln', 'mean', 'dist', and 'randomize' again after filtering for median batch insert size
grouped_data_dynamic = filtered_data_dynamic.groupby(['ln', 'mean', 'dist', 'randomize'])

# Create a figure for the dynamic time plot (for dynamic_gen_time vs num_threads)
plt.figure(figsize=(10, 6))

unique_BI_size = 0;

# Loop through each group
for (ln, mean, dist, randomize), group in grouped_data_dynamic:
    # Check if the configuration has multiple num_threads
    if len(group['num_threads'].unique()) > 1:
        # Filter the group to use dynamic_gen_time and plot against num_threads
        group = group.sort_values('num_threads')
        
        # Generate a unique color for each line
        color = np.random.rand(3,)  # Random RGB color

        unique_BI_size = np.median(group['graph_size']/group['mean'])
        
        # Plot dynamic_gen_time against num_threads
        plt.plot(group['num_threads'], group['dynamic_gen_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}', color=color)

# Customize the dynamic plot
plt.title(f'Dynamic Generation Time vs Number of Threads (batch isnert size = {unique_BI_size})')
plt.xlabel('Number of Threads')
plt.ylabel('Dynamic Generation Time (s)')
plt.grid(True)

# Adjust the legend to be outside, to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Configurations", fontsize=8)

# Save the dynamic time figure
plt.savefig('graphs/dynamic_threads_time.png', bbox_inches='tight')

# Show the dynamic plot
plt.show()
