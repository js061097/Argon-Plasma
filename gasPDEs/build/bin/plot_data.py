
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plot_data.py <X-axis-column> <Y-axis-column>")
    sys.exit(1)

# Get column names for X and Y axes from command-line arguments
x_col = sys.argv[1]
y_col = sys.argv[2]

# Read the CSV file
try:
    df = pd.read_csv('output.csv')
except FileNotFoundError:
    print("Error: 'output.csv' not found. Please run your C++ program to generate the CSV.")
    sys.exit(1)
print(df['t'][0]==df['t'][1])
# Check if the provided columns exist in the CSV
if x_col not in df.columns or y_col not in df.columns:
    print(f"Error: Columns '{x_col}' or '{y_col}' not found in the CSV file.")
    print(f"Available columns: {', '.join(df.columns)}")
    sys.exit(1)

# Plot the data
plt.figure(figsize=(8, 5))
plt.plot(df[x_col], df[y_col], linestyle='-')
plt.xlabel(x_col)
plt.ylabel(y_col)
plt.title(f'{y_col} vs {x_col}')
plt.grid(True)
plt.show()
        