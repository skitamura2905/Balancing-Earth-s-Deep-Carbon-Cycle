import numpy as np
import glob
import os
import matplotlib.pyplot as plt

a = 6371.0e3
b = 3486.0e3
Mman = 4.06e24
Vman = (4 * np.pi / 3) * (a**3 - b**3)  # マントルの体積 (m^3)
Dtra = 660e3
f_man = ((a - Dtra) ** 3 - b**3) / (a**3 - b**3)
Aman = (4 * np.pi) * a**2
g = 9.8

def load_additional_data(file_path="carbon_0.5.dat"):
    """carbon_0.5.dat の 1列目(時間)と 4列目(Vman で正規化)を読み込む"""
    if os.path.exists(file_path):
        data = np.loadtxt(file_path)
        return data[:, 0], data[:, 1], data[:, 2] * 3.15e9
    return None, None, None, None, None

# ファイル読み込み
time, *median_cols = np.loadtxt("thermal_median.dat", unpack=True)
median_cols = np.array(median_cols).T
ci50 = np.load("ci50_thermal.npy")  # (2, N, num_vars)
ci95 = np.load("ci95_thermal.npy")  # (2, N, num_vars)

extra_time, tsurf, uplate = load_additional_data()

# 描画対象インデックス（例：total（マントル合計）は -3 番目）
DATA_SETS = {
    "Tman": ("Mantle Temperature", "mantle_tem.png", "mantle_tem.ps", 0),
    "Tsurf": ("Surface Temperature", "surface_tem.png", "surface_tem.ps", 5),
    "Uplate": ("Plate Velocity", "velocity.png", "velocity.ps", 15),
    "Sp": ("Plate Flux", "slab_volume.png", "slab_volume.ps", 16),
    "Splume": ("Plume Flux", "plume_volume.png", "plume_volume.ps", 17)
}

    
for key, (title, filename, psname, column_index) in DATA_SETS.items():
    fig, ax1 = plt.subplots(figsize=(22, 17))
    #plt.figure(figsize=(9, 7))
    ax1.plot(time, median_cols[:, column_index], label="Median", color="blue")  # 平均値
    ax1.fill_between(time, ci50[0, :, column_index], ci50[1, :, column_index], color="blue", alpha=0.3, label="50% CI")
    ax1.fill_between(time, ci95[0, :, column_index], ci95[1, :, column_index], color="gray", alpha=0.2, label="95% CI")

    if extra_time is not None and tsurf is not None and key in ["Tsurf"]:
        plt.plot(extra_time, tsurf, linestyle="dashed", color="red", label="Foley (2015)")
    if extra_time is not None and uplate is not None and key in ["Uplate"]:
        plt.plot(extra_time, uplate, linestyle="dashed", color="red", label="Foley (2015)")

    xlabel_mapping = {
        "Tman": r'$\mathrm{T_{m}}$',
        "Tsurf": r'$\mathrm{T_{surf}}$',
        "Uplate": r'$\mathrm{u_{plate}}$',
        "Sp": r'$\mathrm{S_{p}}$',
        "Splume": r'$\mathrm{S_{plume}}$'
    }
    xlabel_name = xlabel_mapping.get(key, key)

    if key == "Uplate":
        plt.yscale("log")
        plt.xlabel("Time (Gyr)",fontsize=56)
        plt.ylabel(f"{xlabel_name} (cm/yr)",fontsize=56)
    elif key == "Tsurf":
        plt.xlabel("Time (Gyr)",fontsize=56)
        plt.ylabel(f"{xlabel_name} (K)",fontsize=56)
        plt.axhline(273.15, linestyle="dashed", color="black", label="273.15K")
        #plt.axhline(264.8, linestyle="dotted", color="black", label="264.8K")
        #plt.ylim(273,350)
    elif key == "Tman":
        plt.xlabel("Time (Gyr)",fontsize=56)
        plt.ylabel(f"{xlabel_name} (K)",fontsize=56)
    else:
        plt.yscale("log")
        plt.xlabel("Time (Gyr)",fontsize=56)
        plt.ylabel(f"{xlabel_name} (m³/s)",fontsize=56)
        plt.ylim(5e2,7e4)

    plt.xscale("log")
    plt.tick_params(axis='both', which='major', labelsize=40)
    plt.legend(fontsize=38)
    plt.title(title,fontsize=43)
    plt.savefig(filename, format='png', transparent=True, facecolor='white')
    plt.savefig(psname, format='ps')