# NEPTUNE Graph-Based Two-Sample Testing

This repository provides the core implementation for the **NEPTUNE** test, a graph-based two-sample test using Quantified Ipsen-Mikhailov (QIM) distances and Minimum-Remoteness (MR) matrices. It is designed for comparing groups of network data, such as brain connectomes or other node-aligned network structures.

---

## 🔍 Overview

The NEPTUNE test evaluates group differences between two sets of adjacency matrices representing networks. It combines:

- **Laplacian eigenvalue-based global structure** (Ipsen-Mikhailov distance)
- **Local Hamming distance** (edge-wise difference)
- **Minimum-Remoteness transformation** for high-dimensional data

---

## 🚀 Quick Start

### Prerequisites

- R >= 4.1
- R packages:
  - `igraph`
  - `Rcpp`
  - `igraphdata` (for examples)

### Installation

Clone this repository and source the function file:

```r
# Clone and enter directory
# git clone https://github.com/yourusername/neptune-test.git
# setwd("neptune-test")

# Load the function definitions
source("NEPTUNE_main.R")
```

---

## 📦 Main Function

### `NEPTUNE_Test()`

```r
NEPTUNE_Test(WorkSample, N_Sample1, N_Sample2, n_perm = 1000, delta = 0.5, seed = 11)
```

**Arguments:**

- `WorkSample`: List of adjacency matrices (each an N x N symmetric matrix)
- `N_Sample1`: Number of samples in Group 1
- `N_Sample2`: Number of samples in Group 2
- `n_perm`: Number of permutations for empirical p-value (default 1000)
- `delta`: Correction term for p-value stability (default 0.5)
- `seed`: Random seed for reproducibility

**Returns:**

- `p_QIM`: Empirical p-value based on QIM distances
- `p_MR`: Empirical p-value based on MR transformation
- `Sample_QIM`: Pairwise QIM distance matrix
- `Sample_MR`: Pairwise MR distance matrix

**Output:** Plots MST-kNN clustering results based on QIM distances.

---

## 📊 Distance Metrics

- **QIM Distance:** Combines Hamming and Ipsen-Mikhailov distance between Laplacian spectra.
- **MR Distance:** Derived from QIM using rank-based transformation emphasizing separation.

---

## 📂 File Structure

```
├── NEPTUNE_main.R       # Core functions and NEPTUNE_Test()
├── README.md            # You are here
└── (Optional) data/     # Place example datasets here
```

---

## 🧪 Example Usage

```r
# Example with dummy data
# Create 10 matrices for group 1 and 10 for group 2
sample1 <- replicate(10, matrix(sample(0:1, 100, replace=TRUE), 10), simplify = FALSE)
sample2 <- replicate(10, matrix(sample(0:1, 100, replace=TRUE), 10), simplify = FALSE)
WorkSample <- c(sample1, sample2)

result <- NEPTUNE_Test(WorkSample, N_Sample1 = 10, N_Sample2 = 10)
print(result$p_QIM)
print(result$p_MR)
```

---

## 📜 License

MIT License

---

## 📬 Contact

Created by [Your Name] — feel free to reach out for collaborations or questions.

