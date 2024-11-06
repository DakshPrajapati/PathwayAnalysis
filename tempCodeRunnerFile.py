import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy.interpolate import splrep, splev

# Set limit and index for pathways
limit = 1.8

# Load pathway data
temp = pd.read_csv("C:/coop/pythonScripts/pathwayData.csv")

# Label data based on density cluster
temp['density_cluster'] = np.where(temp['mean_value'].abs() > limit, 'High', 'Low')

# Filtering for high-density genes for labeling
label_data = temp[temp['density_cluster'] == 'High'].groupby('Gene_Phos').tail(1)

# Define custom labels and colors for plotting
custom_labels = {1: "Early G1", 2.5: "Late G1", 6: "S Phase", 9: "G2 Phase", 11.5: "M Phase"}
color_map = {"Low": '#D3D3D3', "High": "purple"}

# Set up the plot
plt.figure(figsize=(10, 6))
texts = []

# Plot Low-density (grey) lines
for gene_phos, group in temp[temp['density_cluster'] == "Low"].groupby('Gene_Phos'):
    # Sort data for proper spline fitting
    group = group.sort_values(by='time')
    # Fit a B-spline to data points
    spline_params = splrep(group['time'], group['value'], s=0.5)  # Adjust smoothness as needed
    smooth_x = np.linspace(group['time'].min(), group['time'].max(), 200)
    smooth_y = splev(smooth_x, spline_params)
    # Plot the line in grey
    plt.plot(smooth_x, smooth_y, color=color_map["Low"], lw=1)

# Plot High-density (purple) lines on top
for gene_phos, group in temp[temp['density_cluster'] == "High"].groupby('Gene_Phos'):
    # Sort data for proper spline fitting
    group = group.sort_values(by='time')
    # Fit a B-spline to data points
    spline_params = splrep(group['time'], group['value'], s=0.5)
    smooth_x = np.linspace(group['time'].min(), group['time'].max(), 200)
    smooth_y = splev(smooth_x, spline_params)
    # Plot the line in purple
    plt.plot(smooth_x, smooth_y, color=color_map["High"], lw=1)
    
    # Add label only at the last point for high-density lines
    last_x = smooth_x[-1]
    last_y = smooth_y[-1]
    texts.append(plt.text(last_x, last_y, gene_phos, color="purple"))

# Set labels and titles
plt.title("Mus musculus: PIP3 activates AKT signaling")  
plt.xlabel("Cell Cycle Phase")
plt.ylabel("Standardized Profile")
plt.xticks(list(custom_labels.keys()), list(custom_labels.values()))
plt.xlim(0, max(custom_labels.keys()) + 1)

# Adjust text to prevent overlapping
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='grey', lw=0.5))

# Show the plot
plt.show()

