import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV data into a pandas DataFrame
data = pd.read_csv('output.csv')

# Group the data by 'ln', 'mean', 'dist', and 'randomize'
grouped_data = data.groupby(['ln', 'mean', 'dist', 'randomize'])

# Create a figure for the plot
plt.figure(figsize=(10, 6))

# Loop through each group and plot
for (ln, mean, dist, randomize), group in grouped_data:
    # Plot the static_gen_time against graph_size for each group
    plt.plot(group['graph_size'], group['static_gen_time'], label=f'ln:{ln} mean:{mean} dist:{dist} randomize:{randomize}')

# Customize the plot
plt.title('Static Generation Time vs Graph Size')
plt.xlabel('Graph Size')
plt.ylabel('Static Generation Time')
plt.legend(loc='best')
plt.grid(True)

# Show the plot
plt.show()
