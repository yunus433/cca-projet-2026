import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

BIT_SIZE = [5]
INCLUDE_RAND = False

def plot_benchmark(mode):
  for bit_size in BIT_SIZE:
    file_name = f'{mode}_expo_bit_size_{bit_size}'

    if not os.path.exists(file_name):
      print(f"ERROR: Could not find {file_name}")
      print("Run benchmark.c to generate data first.")
      return
  
    # Read CSV
    df = pd.read_csv(file_name)
    print(f"Loaded {len(df)} data points")

    # Check required columns
    required_cols = {"length","lib","rand","rand_pos","rand_small","rand_small_pos","order","order_pos"}
    if not required_cols.issubset(df.columns):
      print(f"ERROR: CSV must contain columns: {required_cols}")
      return

    print(df["length"].tolist())

    plt.figure(figsize=(8, 5))
    plt.plot(df["length"].tolist(), df["lib"].tolist(), marker=".", label="Library")
    if INCLUDE_RAND:
      plt.plot(df["length"].tolist(), df["rand"].tolist(), marker=".", label="Random Points")
      plt.plot(df["length"].tolist(), df["rand_pos"].tolist(), marker=".", label="Random Positive Points")
      plt.plot(df["length"].tolist(), df["rand_small"].tolist(), marker=".", label="Small Random Points")
      plt.plot(df["length"].tolist(), df["rand_small_pos"].tolist(), marker=".", label="Small Random Positive Points")
    plt.plot(df["length"].tolist(), df["order"].tolist(), marker=".", label="Ordered Points")
    plt.plot(df["length"].tolist(), df["order_pos"].tolist(), marker=".", label="Ordered Positive Points")
    plt.legend()
    plt.xlabel("Length of the Polynomial")
    plt.ylabel("Execution Time")
    plt.title(f'For {mode} Execution Time vs Length of the Polynomial for Exponent Bitsize {bit_size}')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
  if len(sys.argv) <= 1:
    print("ERROR: Choose mode fmpz or fmpq")
  else:
    mode = sys.argv[1]

    if mode == "fmpz":
      plot_benchmark(mode)
    elif mode == "fmpq":
      plot_benchmark(mode)
    else:
      print("ERROR: Choose mode fmpz or fmpq")