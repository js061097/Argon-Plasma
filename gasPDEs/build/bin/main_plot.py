
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

# Check if the provided columns exist in the CSV
if x_col not in df.columns or y_col not in df.columns:
    print(f"Error: Columns '{x_col}' or '{y_col}' not found in the CSV file.")
    print(f"Available columns: {', '.join(df.columns)}")
    sys.exit(1)
    
time_period = 0.00002

df_filtered = df[df['t']<=time_period]
n_periods = int(df[x_col].iloc[-1]//time_period)
print(df[[x_col,y_col]])
period_vals = []
for i in range(n_periods):
	start = i * time_period
	end = (i+1) * time_period
	period_vals.append(df[(df[x_col]>=start) & (df[x_col]<end)])
print(df[(df[x_col]>=time_period) & (df[x_col]<2*time_period)])

# Plot the data
plt.figure(figsize=(8, 5))
#plt.plot(df[x_col], df[y_col], marker='o', linestyle='-')
plt.plot(df_filtered[x_col], df_filtered[y_col], linestyle='-')
plt.xlabel(x_col)
plt.ylabel(y_col)
plt.title(f'{y_col} vs {x_col}')
plt.grid(True)
plt.show()
        