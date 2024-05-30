import matplotlib.pyplot as plt
import csv
import sys

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

# Main function to plot the data
def plot_data(input_file, plot_title, log_x=False, log_y=False):
    # Read and sort data from CSV file
    num_vertices, time = read_and_sort_data(input_file)
    
    # Plot sorted data
    plt.plot(num_vertices, time, marker='o', linestyle='-')
    
    # Extract the last word from the plot title for the x-axis label
    x_label = plot_title.split()[-1]
    plt.xlabel(x_label)
    
    plt.ylabel('Time (seconds)')
    plt.title(plot_title)
    plt.grid(True)

    # Set x-axis to logarithmic scale if log_x is True
    if log_x:
        plt.xscale('log')

    # Set y-axis to logarithmic scale if log_y is True
    if log_y:
        plt.yscale('log')

    # Save plot as an image file
    output_file = f'{plot_title}.png'
    plt.savefig(output_file)

    # Show plot
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) not in [3, 4, 5]:
        print("Usage: python script.py <input_file> <plot_title> [--log-x] [--log-y]")
    else:
        input_file = sys.argv[1]
        plot_title = sys.argv[2]
        log_x = '--log-x' in sys.argv
        log_y = '--log-y' in sys.argv
        plot_data(input_file, plot_title, log_x, log_y)
