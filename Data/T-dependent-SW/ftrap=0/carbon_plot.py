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
        return data[:, 0], data[:, 3] / Vman, data[:, 5]*Aman*1e5/g, data[:, 6], data[:, 7], data[:, 8]
    return None, None, None, None, None, None

# ファイル読み込み
time, *median_cols = np.loadtxt("carbon_median.dat", unpack=True)
median_cols = np.array(median_cols).T
ci50 = np.load("ci50_carbon.npy")  # (2, N, num_vars)
ci95 = np.load("ci95_carbon.npy")  # (2, N, num_vars)

extra_time, mantle, atm, ocean, sed, bas = load_additional_data()

# 描画対象インデックス（例：total（マントル合計）は -3 番目）
DATA_SETS = {
    "UM C": ("Carbon in Upper Mantle", "upper_mantle.pdf", "upper_mantle.ps", 0),
    "LM C": ("Carbon in Lower Mantle", "lower_mantle.pdf", "lower_mantle.ps", 1),
    "ATM C": ("Carbon in Atmosphere", "atmosphere.pdf", "atmosphere.ps", 2),
    "OC C": ("Carbon in Ocean", "ocean.pdf", "ocean.ps", 3),
    "Sed C": ("Carbon in Ocean Sediment", "sediment.pdf", "sediment.ps", 4),
    "Bas C": ("Carbon in Ocean Basalt", "basalt.pdf", "basalt.ps", 5),
    "Total C": ("Carbon in Whole mantle", "whole_mantle.pdf", "whole_mantle.ps", -3),
    "FD": ("Mid-Ocean Ridge Flux", "mid-ocean.pdf", "mid-ocean.ps", 9),
    "FM": ("Arc Flux", "metamorphic.pdf", "metamorphic.ps", 8),
    "Fplume": ("Plume Flux", "plume.pdf", "plume.ps", 10),
    "Ingas_UM": ("Ingassing Flux to Upper Mantle", "ingas_um.pdf", "ingas_um.ps", 6),
    "Ingas_LM": ("Ingassing Flux to Lower Mantle", "ingas_lm.pdf", "ingas_lm.ps", 7),
    "Inter": ("Lower to Upper Flux", "inter.pdf", "inter.ps", 11),
    "Weathering_c": ("Continental Weathering", "weather_c.pdf", "weather_c.ps", 13),
    "Weathering_s": ("Seafloor Weathering", "weather_s.pdf", "weather_s.ps", 14),
    "Melt_p": ("Plume Residual", "melt_p.pdf", "melt_p.ps", 15),
    "Subd": ("Subduction Flux", "subd.pdf", "subd.ps", -2),
    "Total": ("Total Volcanic Outgassing", "total.pdf", "total.ps", -1)
}

    
for key, (title, filename, psname, column_index) in DATA_SETS.items():
    fig, ax1 = plt.subplots(figsize=(16, 10))
    #plt.figure(figsize=(9, 7))
    ax1.plot(time, median_cols[:, column_index], label="Median", color="blue")  # 平均値
    ax1.fill_between(time, ci50[0, :, column_index], ci50[1, :, column_index], color="blue", alpha=0.3, label="50% CI")
    ax1.fill_between(time, ci95[0, :, column_index], ci95[1, :, column_index], color="gray", alpha=0.2, label="95% CI")

    ax1.set_xscale("log")
    ax1.set_yscale("log")        

    # 追加データのプロット
    if extra_time is not None and mantle is not None:
        if key == "UM C":
            ax2 = ax1.twinx()
            ax2.plot(extra_time, mantle, linestyle="dashed", color="red", label="Foley (2015)")
            ax2.set_ylabel(f"{r'$\mathrm{C_{um}}$'} (kg/m³)", fontsize=35)
            y_min, y_max = ax1.get_ylim()
            ax2.set_ylim(y_min / (Vman * (1 - f_man)), y_max / (Vman * (1 - f_man)))
            ax2.set_yscale("log")
            ax2.tick_params(axis='both', which='major', labelsize=27)
        elif key == "LM C":
            ax2 = ax1.twinx()
            ax2.plot(extra_time, mantle, linestyle="dashed", color="red", label="Foley (2015)")
            ax2.set_ylabel(f"{r'$\mathrm{C_{lm}}$'} (kg/m³)", fontsize=35)
            y_min, y_max = ax1.get_ylim()
            ax2.set_ylim(y_min / (Vman * (f_man)), y_max / (Vman * (f_man)))
            ax2.set_yscale("log")
            ax2.tick_params(axis='both', which='major', labelsize=27)

    if extra_time is not None and atm is not None and key in ["ATM C"]:
        ax1.plot(extra_time, atm, linestyle="dashed", color="red", label="Foley (2015)")
        ax2 = ax1.twinx()
        ax2.set_ylabel(f"{r'$\mathrm{P_{atm}}$'} (bar)", fontsize=35)
        y_min, y_max = ax1.get_ylim()
        ax2.set_ylim(y_min * g/(Aman*1e5), y_max * g/(Aman*1e5))
        ax2.set_yscale("log")    
        ax2.tick_params(axis='both', which='major', labelsize=27)
    if extra_time is not None and ocean is not None and key in ["OC C"]:
        ax1.plot(extra_time, ocean, linestyle="dashed", color="red", label="Foley (2015)")
    if extra_time is not None and sed is not None and key in ["Sed C"]:
        ax1.plot(extra_time, sed, linestyle="dashed", color="red", label="Foley (2015)")
    if extra_time is not None and bas is not None and key in ["Bas C"]:
        ax1.plot(extra_time, bas, linestyle="dashed", color="red", label="Foley (2015)")
        
    # 追加データのプロット
    if key == "UM C" and ax2 is not None:
        target_time = 4.6  # 4.6 Gyr
        target_value = 72e-6 * Mman / Vman  # ppm
        error = 19e-6 * Mman / Vman
        ax2.errorbar(target_time, target_value, yerr=error, fmt='o', color='orange', markersize=8, label='Saal et al. (2002)')
    elif key == "LM C" and ax2 is not None:
        target_time = 4.6  # 4.6 Gyr
        target_value = 263e-6 * Mman / Vman  # ppm
        error_upper = 81e-6 * Mman / Vman
        error_lower = 62e-6 * Mman / Vman
        ax2.errorbar(target_time, target_value, yerr=[[error_lower], [error_upper]], fmt='o', color='green', markersize=8, label='Anderson and Poland (2017)')

    if key == "Total":
        plt.plot(time, median_cols[:, 9], label="Mid-Ocena Ridge Flux", linestyle="dashed")
        plt.plot(time, median_cols[:, 8], label="Arc Flux", linestyle="dashed")
        plt.plot(time, median_cols[:, 10], label="Plume Flux", linestyle="dashed")

    #ax2.set_yscale("log")
    ax1.tick_params(axis='both', which='major', labelsize=27)
    plt.title(title,fontsize=32)
        
    xlabel_mapping = {
        "UM C": r'$\mathrm{M_{um}}$',
        "LM C": r'$\mathrm{M_{lm}}$',
        "ATM C": r'$\mathrm{M_{atm}}$',
        "OC C": r'$\mathrm{M_{oc}}$',
        "Sed C": r'$\mathrm{M_{sed}}$',
        "Bas C": r'$\mathrm{M_{bas}}$',
        "Total C": r'$\mathrm{M_{m}}$',
        "FD": r'$\mathrm{F_{mor}}$',
        "FM": r'$\mathrm{F_{arc}}$',
        "Fplume": r'$\mathrm{F_{plume}}$',
        "Ingas_UM": r'$\mathrm{f_{slab}F_{ingas}}$',
        "Ingas_LM": r'$\mathrm{(1-f_{slab})F_{ingas}}$',
        "Inter": r'$\mathrm{F_{inter}}$',
        "Weathering_c": r'$\mathrm{F_{sil}}$',
        "Weathering_s": r'$\mathrm{F_{sw}}$',
        "Melt_mor": r'$\mathrm{(1/\epsilon_{erupt,mor}-1)F_{mor}}$',
        "Melt_p": r'$\mathrm{F_{res}}$',
        "Melt_lower": r'$\mathrm{f_{trap}F_{trap}}$',
        "Subd": r'$\mathrm{F_{ingas}}$',
        "Total": r'$\mathrm{F_{total}}$'
    }
    xlabel_name = xlabel_mapping.get(key, key)

    # x軸の単位を変更  
    if key in ["OC C", "Sed C", "Bas C", "ATM C", "UM C", "LM C", "Total C"]:
        ax1.set_xlabel("Time (Gyr)",fontsize=35)
        ax1.set_ylabel(f"{xlabel_name} (kg)",fontsize=35)
    else:
        ax1.set_xlabel("Time (Gyr)",fontsize=35)
        ax1.set_ylabel(f"{xlabel_name} (kg/yr)",fontsize=35)
        ax1.set_ylim(1e9,1e15)
        

    if key in ["UM C", "LM C", "ATM C"]:
        lines, labels = ax1.get_legend_handles_labels()
        if ax2 is not None:
            lines2, labels2 = ax2.get_legend_handles_labels()
            lines += lines2
            labels += labels2
            # 1つの凡例に統合
            ax1.legend(lines, labels, fontsize=20, loc="best")
    else:
        ax1.legend(fontsize=20)

    plt.savefig(filename, format="pdf", transparent=True, facecolor="white")
    plt.savefig(psname, format="ps")
