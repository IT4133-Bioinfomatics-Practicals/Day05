# 🧬 Day 05: Gene Expression Data Analysis & Visualization Guide

📊 A complete guide to analyzing and visualizing differential gene expression (DEG) data with practical Python examples.

---

## 📋 Overview

This notebook analyzes gene expression data comparing **Control** vs **Treatment** conditions across multiple samples. It includes calculations for average expression, identification of regulated genes, and 7 different visualization techniques.

---

## 📁 Dataset Structure

| Column                    | Description                                  |
| ------------------------- | -------------------------------------------- |
| 🧬 Gene                   | Gene name (Gene1, Gene2, etc.)               |
| 🟦 Control1-Control9      | Expression values under control conditions   |
| 🟥 Treatment1-Treatment25 | Expression values under treatment conditions |

---

## 💻 Code Breakdown

### 1️⃣ **Load & Prepare Data**

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("dea_data.csv")
```

- Imports essential libraries for data analysis and visualization
- Loads the CSV file into a pandas DataFrame

### 2️⃣ **Extract Column Names**

```python
control_cols = [i for i in data.columns if "Control" in i]
treatment_cols = [i for i in data.columns if "Treatment" in i]
```

- Automatically identifies control and treatment columns
- Creates lists for easy reference in calculations

### 3️⃣ **Calculate Average Expression (Control)**

```python
for gene in data["Gene"].head(100):
    gene_data = data[data["Gene"] == gene]
    control_mean = gene_data[control_cols].mean(axis=1).values[0]
    print(f"Gene: {gene}, Average: {control_mean:.2f}")
```

- Loops through first 100 genes
- Calculates mean expression for each gene under control conditions

### 4️⃣ **Calculate Average Expression (Treatment)**

```python
for gene in data["Gene"].head(100):
    gene_data = data[data["Gene"] == gene]
    treatment_mean = gene_data[treatment_cols].mean(axis=1).values[0]
    print(f"Gene: {gene}, Average: {treatment_mean:.2f}")
```

- Similar to above but for treatment conditions

### 5️⃣ **Find Highest Expressed Gene (Treatment)**

```python
data["treatment_avg"] = data[treatment_cols].mean(axis=1)
highest_gene = data.loc[data["treatment_avg"].idxmax()]
print(f"Gene: {highest_gene['Gene']}, Avg: {highest_gene['treatment_avg']:.2f}")
```

- Creates a new column with treatment averages
- Identifies gene with highest average expression

### 6️⃣ **Calculate Expression Difference**

```python
data["control_avg"] = data[control_cols].mean(axis=1)
data["expression_diff"] = data["treatment_avg"] - data["control_avg"]
print(data[["Gene", "control_avg", "treatment_avg", "expression_diff"]])
```

- Computes average for control conditions
- Calculates the difference: Treatment - Control

### 7️⃣ **Bar Chart Visualization** 📊

```python
gene = data["Gene"].iloc[0]  # First gene
gene_data = data[data["Gene"] == gene]
control_mean = gene_data[control_cols].mean(axis=1).values[0]
treatment_mean = gene_data[treatment_cols].mean(axis=1).values[0]

plt.bar("Control", control_mean, color='blue', label='Control')
plt.bar("Treatment", treatment_mean, color='orange', label='Treatment')
plt.title(f"Bar Chart: {gene} Average Expression")
plt.ylabel("Average Expression Level")
plt.legend()
plt.show()
```

- **Purpose:** Compare control vs treatment for a single gene
- **Best for:** Quick comparison of two conditions

### 8️⃣ **Heatmap Visualization** 🔥

```python
heatmap_data = data.set_index("Gene").head(100)
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, cmap="coolwarm", annot=True)
plt.title("Heatmap of Gene Expression Levels")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.show()
```

- **Purpose:** View expression patterns across all genes and samples
- **Best for:** Identifying clusters and patterns

### 9️⃣ **Identify Upregulated Genes** ⬆️

```python
upregulated_gene = data[data["treatment_avg"] > data["control_avg"]]
print(upregulated_gene[["Gene", "control_avg", "treatment_avg", "expression_diff"]])
```

- **Purpose:** Find genes with higher expression in treatment
- **Filter:** `treatment_avg > control_avg`

### 🔟 **Identify Downregulated Genes** ⬇️

```python
downregulated_gene = data[data["treatment_avg"] < data["control_avg"]]
print(downregulated_gene[["Gene", "control_avg", "treatment_avg", "expression_diff"]])
```

- **Purpose:** Find genes with lower expression in treatment
- **Filter:** `treatment_avg < control_avg`

### 1️⃣1️⃣ **Box Plot for First 5 Genes** 📦

```python
first_5_genes = data.index[:5]
treatment_data_5 = data.loc[first_5_genes, treatment_cols]
plt.figure(figsize=(10, 6))
treatment_data_5.T.boxplot(rot=45)
plt.title('Treatment Expression Levels for First 5 Genes')
plt.ylabel('Expression Level')
plt.xlabel('Genes')
plt.tight_layout()
plt.show()
```

- **Purpose:** Show distribution and outliers in treatment data
- **Best for:** Identifying variability and outliers

### 1️⃣2️⃣ **Overall Mean Expression** 📈

```python
overall_control_mean = data[control_cols].mean().mean()
overall_treatment_mean = data[treatment_cols].mean().mean()
print(f"Control: {overall_control_mean}, Treatment: {overall_treatment_mean}")
```

- **Purpose:** Compare global expression levels
- **Calculation:** Mean of all control + mean of all treatment

### 1️⃣3️⃣ **Find Most Variable Gene** 🎲

```python
all_samples = control_cols + treatment_cols
data["Variation"] = data[all_samples].std(axis=1)
most_variable_gene = data.loc[data["Variation"].idxmax()]
print(f"Gene: {most_variable_gene['Gene']}, Variation: {most_variable_gene['Variation']:.2f}")
```

- **Purpose:** Identify gene with highest variability
- **Metric:** Standard deviation across all samples

### 1️⃣4️⃣ **Scatter Plot: Upregulated vs Downregulated** 🔴🔵

```python
plt.figure(figsize=(8, 6))
colors = {'upregulate':'red', 'downregulate':'blue'}
summary_df = pd.DataFrame({
    'control_avg': data[control_cols].mean(axis=1),
    'treatment_avg': data[treatment_cols].mean(axis=1)
})
summary_df['Regulation'] = summary_df.apply(
    lambda row: 'upregulate' if row['treatment_avg'] > row['control_avg'] else 'downregulate',
    axis=1
)

for reg in ['upregulate', 'downregulate']:
    subset = summary_df[summary_df['Regulation'] == reg]
    plt.scatter(subset['control_avg'], subset['treatment_avg'],
                c=colors[reg], label=reg, alpha=0.7)

plt.plot([30, summary_df[['control_avg', 'treatment_avg']].max().max()],
         [30, summary_df[['control_avg', 'treatment_avg']].max().max()], 'k--', label='y=x')
plt.xlabel("Control Average Expression")
plt.ylabel("Treatment Average Expression")
plt.title("Upregulate vs Downregulate Genes")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

- **Purpose:** Visualize regulation status of all genes
- **Red points:** Upregulated (treatment > control)
- **Blue points:** Downregulated (treatment < control)
- **Diagonal line:** y=x (no change)

---

## 🔑 Key Concepts

| Term                   | Definition                                                                |
| ---------------------- | ------------------------------------------------------------------------- |
| ⬆️ **Upregulated**     | Gene expression increases in treatment (treatment_avg > control_avg)      |
| ⬇️ **Downregulated**   | Gene expression decreases in treatment (treatment_avg < control_avg)      |
| 🧬 **DEG**             | Differentially Expressed Gene - genes with significant expression changes |
| 📊 **Mean Expression** | Average value across all samples for a condition                          |
| 📈 **Fold Change**     | Ratio of treatment expression to control expression                       |

---

## 💡 Quick Tips

✅ **Always run cells in order** - Each cell depends on previous calculations
✅ **Use `.head(n)`** - Test with small datasets first (e.g., `.head(100)`)
✅ **Check column names** - Use `data.columns` if unsure about naming
✅ **Color coding** - Use consistent colors (🔵=control, 🔴=treatment)
✅ **Add labels** - Always include title, xlabel, ylabel in plots

---

## ⚠️ Common Errors & Fixes

| ❌ Error                               | 🔍 Cause                  | ✅ Fix                                                 |
| -------------------------------------- | ------------------------- | ------------------------------------------------------ |
| `KeyError: 'control_avg'`              | Column doesn't exist      | Run cells in order; create column first                |
| `IndexError: index 0 is out of bounds` | Gene not found            | Use `data["Gene"].iloc[0]` instead of hardcoded names  |
| `FileNotFoundError`                    | CSV file not in directory | Check file path or use `pd.read_csv("./dea_data.csv")` |

---

## 🐼 Useful Pandas Functions

```python
data.head(n)           # First n rows
data.tail(n)           # Last n rows
data[col].mean()       # Average value
data[col].std()        # Standard deviation
data[col].min()        # Minimum value
data[col].max()        # Maximum value
data.loc[condition]    # Filter rows by condition
data.set_index(col)    # Set column as index
```

---

## 🎨 Matplotlib & Seaborn Cheat Sheet

```python
plt.figure(figsize=(width, height))    # Set figure size
plt.bar(x, y)                          # Bar chart
plt.scatter(x, y)                      # Scatter plot
plt.plot(x, y)                         # Line plot
plt.title("Title")                     # Add title
plt.xlabel("X Label")                  # X-axis label
plt.ylabel("Y Label")                  # Y-axis label
plt.legend()                           # Show legend
plt.show()                             # Display plot
sns.heatmap(data)                      # Heatmap
```

---

📅 **Created:** November 2025  
🔬 **Topic:** Bioinformatics - Gene Expression Analysis  
🎓 **Course:** IT 4133 - Bioinfomatics and Computational Biology
