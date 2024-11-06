# Cell Cycle Pathway Visualization

## Overview
This Python script generates a visualization of gene phosphorylation pathways across different cell cycle phases. It distinguishes between high and low-density pathways using smooth spline interpolation and creates an intuitive plot with labeled significant pathways.

## Features
- Visualizes gene phosphorylation patterns across the cell cycle
- Distinguishes between high and low-density pathways using a density threshold
- Implements smooth spline interpolation for better curve visualization
- Automatically labels high-density pathways with non-overlapping text
- Color-codes pathways (purple for high-density, grey for low-density)
- Marks distinct cell cycle phases (Early G1, Late G1, S Phase, G2 Phase, M Phase)

## Prerequisites
The script requires the following Python packages:
```python
pandas
numpy
seaborn
matplotlib
adjustText
scipy