import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import string # アルファベット生成用

# --- Constants ---
a = 6371.0e3
b = 3486.0e3
Mman = 4.06e24
Vman = (4 * np.pi / 3) * (a**3 - b**3)  # マントルの体積 (m^3)
Dtra = 660e3
f_man = ((a - Dtra) ** 3 - b**3) / (a**3 - b**3)
Aman = (4 * np.pi) * a**2
g = 9.8

# アルファベット生成用のイテレータ
alphabet_iterator = iter(string.ascii_lowercase)

# --- Data Loading Functions ---

def load_carbon_additional_data(file_path="carbon_0.5.dat"):
    """
    Loads carbon-related additional data from carbon_0.5.dat.
    Returns: time, mantle (normalized by Vman), atm (normalized), ocean, plate
    """
    if os.path.exists(file_path):
        data = np.loadtxt(file_path)
        # Columns: 0:time, 1:tsurf, 2:uplate, 3:mantle, 4:atm (mass), 5:atm (pres), 6:ocean, 7:plate
        return data[:, 0], data[:, 3], data[:, 5]*Aman*1e5/g, data[:, 6], data[:, 7], data[:, 8], data[:, 9], data[:, 10], data[:, 11], data[:, 12]
    return None, None, None, None, None, None, None, None, None

def load_thermal_additional_data(file_path="carbon_0.5.dat"):
    """
    Loads thermal-related additional data from carbon_0.5.dat (assuming it contains Tsurf and Uplate).
    Returns: time, tsurf, uplate
    """
    if os.path.exists(file_path):
        data = np.loadtxt(file_path)
        # Columns: 0:time, 1:tsurf, 2:uplate, 3:mantle, 4:atm (mass), 5:atm (pres), 6:ocean, 7:plate
        return data[:, 0], data[:, 1], data[:, 2] * 3.15e9 # 3.15e9 is likely for unit conversion (e.g., cm/s to cm/yr)
    return None, None, None


# --- Plotting Configuration ---

# Data set definitions for carbon-related plots
DATA_SETS_CARBON = {
    "M C": ("Carbon in Mantle", "mantle.png", "mantle.ps", 0),
    "ATM C": ("Carbon in Atmosphere", "atmosphere.png", "atmosphere.ps", 2),
    "OC C": ("Carbon in Ocean", "ocean.png", "ocean.ps", 3),
    "P C": ("Carbon in Ocean Crust", "plate.png", "plate.ps", 4),
    "FD": ("Mid-Ocean Ridge Flux", "mid-ocean.png", "mid-ocean.ps", 7),
    "FM": ("Arc Flux", "metamorphic.png", "metamorphic.ps", 6),
    "Ingas": ("Ingassing Flux", "ingas_um.png", "ingas_um.ps", 5),
    "Weathering_c": ("Continental Weathering", "weather_c.png", "weather_c.ps", 8),
    "Weathering_s": ("Seafloor Weathering", "weather_s.png", "weather_s.ps", 9)
}

# Data set definitions for thermal-related plots
DATA_SETS_THERMAL = {
    "Tman": ("Mantle Temperature", "mantle_tem.png", "mantle_tem.pdf", 0),
    "Tsurf": ("Surface Temperature", "surface_tem.png", "surface_tem.pdf", 4),
    "Uplate": ("Plate Velocity", "velocity.png", "velocity.pdf", 12),
    "Sp": ("Plate Flux", "slab_volume.png", "slab_volume.pdf", 13)
}

# LaTeX label mappings for Y-axis
xlabel_mapping = {
    "M C": r'$\mathrm{M_{m}}$',
    "ATM C": r'$\mathrm{M_{atm}}$', "OC C": r'$\mathrm{M_{oc}}$',
    "P C": r'$\mathrm{M_{plate}}$',
    "FD": r'$\mathrm{F_{mor}}$', "FM": r'$\mathrm{F_{arc}}$',
    "Ingas": r'$\mathrm{F_{ingas}}$',
    "Weathering_c": r'$\mathrm{F_{sil}}$', "Weathering_s": r'$\mathrm{F_{sw}}$',
    "Tman": r'$\mathrm{T_{m}}$', "Tsurf": r'$\mathrm{T_{surf}}$',
    "Uplate": r'$\mathrm{u_{plate}}$', "Sp": r'$\mathrm{S_{p}}$'
}

# --- Load Data (once for all plots) ---
time_carbon, *median_cols_carbon_list = np.loadtxt("carbon_median.dat", unpack=True)
median_cols_carbon = np.array(median_cols_carbon_list).T
ci50_carbon = np.load("ci50_carbon.npy")  # (2, N, num_vars)
ci95_carbon = np.load("ci95_carbon.npy")  # (2, N, num_vars)
extra_time_carbon, mantle_extra, atm_extra, ocean_extra, plate_extra, ingas, meta, degas, fcw, fsw = load_carbon_additional_data()

time_thermal, *median_cols_thermal_list = np.loadtxt("thermal_median.dat", unpack=True)
median_cols_thermal = np.array(median_cols_thermal_list).T
ci50_thermal = np.load("ci50_thermal.npy")  # (2, N, num_vars)
ci95_thermal = np.load("ci95_thermal.npy")  # (2, N, num_vars)
extra_time_thermal, tsurf_extra, uplate_extra = load_thermal_additional_data()


# --- Plotting Carbon Reservoir and Surface Temperature Data in 2 Columns ---
combined_keys = ["M C", "P C", "ATM C", "OC C", "Tsurf"]
num_plots = len(combined_keys)
nrows = (num_plots + 1) // 2
ncols = 2

fig_width_per_col = 35
fig_height_per_row = 20
fig_combined, axes_combined = plt.subplots(nrows, ncols, figsize=(fig_width_per_col * ncols, fig_height_per_row * nrows))
axes_combined = axes_combined.flatten()

# === ここを修正しました: subplots_adjust で上下の余白をさらに減らす ===
# top: 上端までのサブプロットの終了位置 (0.0 - 1.0)
# bottom: 下端からのサブプロットの開始位置 (0.0 - 1.0)
plt.subplots_adjust(left=0.07, right=0.93, top=0.96, bottom=0.06, wspace=0.4, hspace=0.3)

# アルファベット生成用のイテレータ
alphabet_iterator = iter(string.ascii_lowercase)

print("Generating Combined Carbon and Surface Temperature Multi-plot with minimal vertical margins...")
for i, key in enumerate(combined_keys):
    ax1 = axes_combined[i]
    
    # Determine which dataset to use (carbon or thermal)
    if key in DATA_SETS_CARBON:
        data_set = DATA_SETS_CARBON
        current_time = time_carbon
        median_cols = median_cols_carbon
        ci50 = ci50_carbon
        ci95 = ci95_carbon
        current_extra_time = extra_time_carbon
    elif key in DATA_SETS_THERMAL:
        data_set = DATA_SETS_THERMAL
        current_time = time_thermal
        median_cols = median_cols_thermal
        ci50 = ci50_thermal
        ci95 = ci95_thermal
        current_extra_time = extra_time_thermal
    else:
        print(f"Warning: Key '{key}' not found in any DATA_SETS. Skipping.")
        continue

    title, _, _, column_index = data_set[key]

    # Plot median and confidence intervals
    ax1.plot(current_time, median_cols[:, column_index], label="Median", color="blue", linewidth=6)
    ax1.fill_between(current_time, ci50[0, :, column_index], ci50[1, :, column_index], color="blue", alpha=0.3, label="50% CI")
    ax1.fill_between(current_time, ci95[0, :, column_index], ci95[1, :, column_index], color="gray", alpha=0.2, label="95% CI")

    ax1.set_xscale("log")
    
    # Specific y-axis scaling and additional data for each key
    ax2 = None
    if key in ["M C", "P C", "ATM C", "OC C"]:
        ax1.set_yscale("log")
        if current_extra_time is not None:
            if key == "M C" and mantle_extra is not None:
                ax1.plot(current_extra_time, mantle_extra, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)

            elif key == "ATM C" and atm_extra is not None:
                ax1.plot(current_extra_time, atm_extra, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
                ax2 = ax1.twinx()
                ax2.set_ylabel(r'$\mathrm{P_{atm}}$ (bar)', fontsize=100)
                y_min, y_max = ax1.get_ylim()
                ax2.set_ylim(y_min * g/(Aman*1e5), y_max * g/(Aman*1e5))
                ax2.set_yscale("log")
                ax2.tick_params(axis='both', which='major', labelsize=80)

            elif key == "OC C" and ocean_extra is not None:
                ax1.plot(current_extra_time, ocean_extra, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
            elif key == "P C" and plate_extra is not None:
                ax1.plot(current_extra_time, plate_extra, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
        
        ax1.set_ylabel(f"{xlabel_mapping.get(key, key)} (kg)", fontsize=100)

    elif key == "Tsurf":
        ax1.set_ylabel(f"{xlabel_mapping.get(key, key)} (K)", fontsize=100)
        ax1.axhline(273.15, linestyle="dashed", color="black", label="273.15K")
        if current_extra_time is not None and tsurf_extra is not None:
            ax1.plot(current_extra_time, tsurf_extra, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)

    # Set common x-axis label and tick parameters
    ax1.set_xlabel("Time (Gyr)", fontsize=100)
    ax1.tick_params(axis='x', which='major', labelsize=80, pad=20)
    ax1.tick_params(axis='y', which='major', labelsize=80)
    ax1.set_title(title, fontsize=90)

    # Handle legends for twinx plots
    if ax2 is not None:
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=50, loc="best")
    else:
        ax1.legend(fontsize=50, loc="best")

    # === ここにアルファベットの追加 ===
    try:
        label_char = next(alphabet_iterator)
        ax1.text(-0.17, 1.1, f'{label_char}', transform=ax1.transAxes,
                fontsize=100, va='top', ha='left', weight='bold') # フォントサイズと太さを調整
    except StopIteration:
        # アルファベットが尽きた場合（26個以上サブプロットがある場合）
        pass

# Hide any unused subplots
for i in range(num_plots, nrows * ncols):
    fig_combined.delaxes(axes_combined[i])

# plt.tight_layout() をコメントアウトまたは削除し、subplots_adjust を優先
# plt.tight_layout()

plt.savefig("reservoir.png", format="png", transparent=True, facecolor="white")
plt.savefig("reservoir.pdf", format="pdf", transparent=True, facecolor="white")
plt.close(fig_combined)
print("Combined Carbon and Surface Temperature Multi-plot with minimal vertical margins generated successfully.")


# --- Plotting Carbon Reservoir and Surface Temperature Data in 2 Columns ---
combined_keys = ["FD", "FM", "Ingas", "Weathering_c", "Weathering_s"]
num_plots = len(combined_keys)
nrows = (num_plots + 1) // 2
ncols = 2

fig_width_per_col = 35
fig_height_per_row = 20
fig_combined, axes_combined = plt.subplots(nrows, ncols, figsize=(fig_width_per_col * ncols, fig_height_per_row * nrows))
axes_combined = axes_combined.flatten()

# === ここを修正しました: subplots_adjust で上下の余白をさらに減らす ===
# top: 上端までのサブプロットの終了位置 (0.0 - 1.0)
# bottom: 下端からのサブプロットの開始位置 (0.0 - 1.0)
plt.subplots_adjust(left=0.07, right=0.93, top=0.96, bottom=0.06, wspace=0.4, hspace=0.3)

# アルファベット生成用のイテレータ
alphabet_iterator = iter(string.ascii_lowercase)

print("Generating Combined Carbon and Surface Temperature Multi-plot with minimal vertical margins...")
for i, key in enumerate(combined_keys):
    ax1 = axes_combined[i]
    
    # Determine which dataset to use (carbon or thermal)
    if key in DATA_SETS_CARBON:
        data_set = DATA_SETS_CARBON
        current_time = time_carbon
        median_cols = median_cols_carbon
        ci50 = ci50_carbon
        ci95 = ci95_carbon
        current_extra_time = extra_time_carbon
    elif key in DATA_SETS_THERMAL:
        data_set = DATA_SETS_THERMAL
        current_time = time_thermal
        median_cols = median_cols_thermal
        ci50 = ci50_thermal
        ci95 = ci95_thermal
        current_extra_time = extra_time_thermal
    else:
        print(f"Warning: Key '{key}' not found in any DATA_SETS. Skipping.")
        continue

    title, _, _, column_index = data_set[key]

    # Plot median and confidence intervals
    ax1.plot(current_time, median_cols[:, column_index], label="Median", color="blue", linewidth=6)
    ax1.fill_between(current_time, ci50[0, :, column_index], ci50[1, :, column_index], color="blue", alpha=0.3, label="50% CI")
    ax1.fill_between(current_time, ci95[0, :, column_index], ci95[1, :, column_index], color="gray", alpha=0.2, label="95% CI")

    if current_extra_time is not None:
        if key == "FD" and mantle_extra is not None:
            ax1.plot(current_extra_time, degas, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
        elif key == "FM" and mantle_extra is not None:
            ax1.plot(current_extra_time, meta, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
        elif key == "Weathering_c" and ocean_extra is not None:
            ax1.plot(current_extra_time, fcw, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
        elif key == "Weathering_s" and ocean_extra is not None:
            ax1.plot(current_extra_time, fsw, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)
        elif key == "Ingas" and ocean_extra is not None:
            ax1.plot(current_extra_time, ingas, linestyle="dashed", color="red", label="Foley (2015)", linewidth=6)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylim(1e9,1e15)
    ax1.set_ylabel(f"{xlabel_mapping.get(key, key)} (kg/yr)", fontsize=100)

    # Set common x-axis label and tick parameters
    ax1.set_xlabel("Time (Gyr)", fontsize=100)
    ax1.tick_params(axis='x', which='major', labelsize=80, pad=20)
    ax1.tick_params(axis='y', which='major', labelsize=80)
    ax1.set_title(title, fontsize=90)

    # Handle legends for twinx plots
    if ax2 is not None:
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=50, loc="best")
    else:
        ax1.legend(fontsize=50, loc="best")

    # === ここにアルファベットの追加 ===
    try:
        label_char = next(alphabet_iterator)
        ax1.text(-0.17, 1.1, f'{label_char}', transform=ax1.transAxes,
                fontsize=100, va='top', ha='left', weight='bold') # フォントサイズと太さを調整
    except StopIteration:
        # アルファベットが尽きた場合（26個以上サブプロットがある場合）
        pass

# Hide any unused subplots
for i in range(num_plots, nrows * ncols):
    fig_combined.delaxes(axes_combined[i])

# plt.tight_layout() をコメントアウトまたは削除し、subplots_adjust を優先
# plt.tight_layout()

plt.savefig("flux.png", format="png", transparent=True, facecolor="white")
plt.savefig("flux.pdf", format="pdf", transparent=True, facecolor="white")
plt.close(fig_combined)
print("Combined Carbon and Surface Temperature Multi-plot with minimal vertical margins generated successfully.")

print("All combined plots generated successfully.")