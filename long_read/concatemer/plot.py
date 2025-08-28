import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

#path
t0_path = "/concatemer/result_concatemer/t0.txt"
t120_path = "/concatemer/result_concatemer/t120.txt"
save_dir = "/concatemer/result_concatemer/"

# read files
t0 = pd.read_csv(t0_path, header=None, names=["length"])
t120 = pd.read_csv(t120_path, header=None, names=["length"])

# function for graphs
def plot_kde(data, title, xlim=None, filename="output.png"):
    mean_val = data["length"].mean()
    plt.figure(figsize=(10, 5))
    sns.kdeplot(data["length"], bw_adjust=0.5, fill=True)
    plt.axvline(mean_val, color='red', linestyle='--', label=f"mean = {mean_val:.1f}")
    plt.title(title)
    plt.xlabel("Softclip Length")
    plt.ylabel("Density")
    if xlim:
        plt.xlim(xlim)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path, dpi=300)
    plt.show()
    print(f" Saved: {save_path}")

# graph1- all 
plot_kde(t0, "T0 - Softclip Lengths (all)", filename="t0_all.png")
plot_kde(t120, "T120 - Softclip Lengths (all)", filename="t120_all.png")

# graph2– between 0 and 5000 
plot_kde(t0[t0["length"] <= 5000], "T0 - Softclip Lengths (0–5000)", xlim=(0, 5000), filename="t0_0_5000.png")
plot_kde(t120[t120["length"] <= 5000], "T120 - Softclip Lengths (0–5000)", xlim=(0, 5000), filename="t120_0_5000.png")

# garph 3 – 200 and  5000
plot_kde(t0[(t0["length"] >= 200) & (t0["length"] <= 5000)], "T0 - Softclip Lengths (200–5000)", xlim=(200, 5000), filename="t0_200_5000.png")
plot_kde(t120[(t120["length"] >= 200) & (t120["length"] <= 5000)], "T120 - Softclip Lengths (200–5000)", xlim=(200, 5000), filename="t120_200_5000.png")
