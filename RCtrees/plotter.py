# All this code is copied from chatGPT
# I will not pretend to own it even if this code is pretty trivial
import matplotlib.pyplot as plt
import csv

# Function to read and sort data from CSV file
def read_and_sort_data(filename):
    data = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row
        for row in reader:
            data.append((int(row[0]), float(row[1])))
    sorted_data = sorted(data, key=lambda x: x[0])  # Sort by num_vertices
    return zip(*sorted_data)

# Read and sort data from CSV file
filename = 'data.csv'
num_vertices, time = read_and_sort_data(filename)

# Plot sorted data
plt.plot(num_vertices, time, marker='o', linestyle='-')
plt.xlabel('Number of Vertices')
plt.ylabel('Time (seconds)')
plt.title('Time vs. Number of Vertices (Sorted)')
plt.grid(True)

# Save plot as an image file
plt.savefig('time_vs_vertices_sorted.png')

# Show plot
plt.show()
