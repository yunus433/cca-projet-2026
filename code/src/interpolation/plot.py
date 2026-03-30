import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

def plot_benchmark_dependence(csv_file='data-dependence.csv'):
  # Check if file exists
  if not os.path.exists(csv_file):
    print(f"ERROR: Could not find {csv_file}")
    print("Run benchmark.c to generate data first.")
    return
  
  # Read CSV
  df = pd.read_csv(csv_file)
  print(f"Loaded {len(df)} data points")

  # Check required columns
  required_cols = {"degree", "n", "time"}
  if not required_cols.issubset(df.columns):
    print(f"ERROR: CSV must contain columns: {required_cols}")
    return

  # Remove invalid rows for logarithm
  df = df[(df["n"] > 0) & (df["time"] > 0)].copy()

  results = []

  # Fit log(time) = log(a_d) + k_d * log(n) for each degree
  for degree, group in df.groupby("degree"):
    if len(group) < 2:
      print(f"INFO: Skipping degree {degree}: not enough points for regression")
      continue

    x = np.log(group["n"].values)
    y = np.log(group["time"].values)

    # Linear regression in log-log space
    k_d, intercept = np.polyfit(x, y, 1)
    a_d = np.exp(intercept)

    results.append({
        "degree": degree,
        "a_d": a_d,
        "k_d": k_d
    })

  print("here")

  results_df = pd.DataFrame(results).sort_values("degree")

  print(results_df)

  # Plot a_d vs degree
  plt.figure(figsize=(8, 5))
  plt.plot(results_df["degree"], results_df["a_d"], marker=".")
  plt.xlabel("Degree")
  plt.ylabel("a_d")
  plt.title("Coefficient a_d vs Degree")
  plt.grid(True)
  plt.tight_layout()
  plt.show()

  # Plot k_d vs degree
  plt.figure(figsize=(8, 5))
  plt.plot(results_df["degree"], results_df["k_d"], marker=".")
  plt.xlabel("Degree")
  plt.ylabel("k_d")
  plt.title("Exponent k_d vs Degree")
  plt.grid(True)
  plt.tight_layout()
  plt.show()

def plot_benchmark_pow_2(csv_file='data-pow-2.csv'):
  # Check if file exists
  if not os.path.exists(csv_file):
    print(f"ERROR: Could not find {csv_file}")
    print("Run benchmark.c to generate data first.")
    return
  
  # Read CSV
  df = pd.read_csv(csv_file)
  print(f"Loaded {len(df)} data points")

  # Check required columns
  required_cols = {"degree", "ratio"}
  if not required_cols.issubset(df.columns):
    print(f"ERROR: CSV must contain columns: {required_cols}")
    return

  # Plot ratio vs degree
  plt.figure(figsize=(8, 5))
  plt.plot(df["degree"], df["ratio"], marker=".")
  plt.xlabel("Degree")
  plt.ylabel("Ratio (pow2 / normal)")
  plt.title("Ratio over Degree")
  plt.grid(True)
  plt.tight_layout()
  plt.show()

if __name__ == '__main__':
  if len(sys.argv) <= 1:
    print("ERROR: Choose mode DEPENDENCE or POW2")
  else:
    mode = sys.argv[1]

    if mode == "DEPENDENCE":
      csv_file = sys.argv[2] if len(sys.argv) > 2 else 'data-dependence.csv'
      plot_benchmark_dependence(csv_file)
    elif mode == "POW2":
      csv_file = sys.argv[2] if len(sys.argv) > 2 else 'data-pow-2.csv'
      plot_benchmark_pow_2(csv_file)
    else:
      print("ERROR: Choose mode DEPENDENCE or POW2")
  