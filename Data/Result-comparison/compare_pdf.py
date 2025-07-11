import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages # PDF出力用
import string # アルファベット生成用

# --- 定数 ---
a = 6371.0e3
b = 3486.0e3
Mman = 4.06e24
Vman = (4 * np.pi / 3) * (a**3 - b**3)  # マントルの体積 (m^3)
Dtra = 660e3
f_man = ((a - Dtra) ** 3 - b**3) / (a**3 - b**3)
Aman = 4 * np.pi * a**2 # Aman の計算が前回間違っていたため修正
g = 9.8

# --- データセット情報と列インデックス ---
# 各キーに対応する: (タイトル, 列インデックス, 単位, ylim_min, ylim_max, y2_label (オプション))
# 列インデックスはcarbon_median_standard.datを基準とする
# 注: melt_lower と Fres はデータファイルによって列インデックスが異なるため、特別扱いが必要
DATA_CONFIG = {
    "UM C":        {"title": "Carbon in Upper Mantle",      "col_idx": 1,  "unit": "(kg)",       "y2_label": r'$\mathrm{C_{um}}$ (kg/m³)', "y2_scale_factor": lambda y_val: y_val / (Vman * (1 - f_man))},
    "LM C":        {"title": "Carbon in Lower Mantle",      "col_idx": 2,  "unit": "(kg)",       "y2_label": r'$\mathrm{C_{lm}}$ (kg/m³)', "y2_scale_factor": lambda y_val: y_val / (Vman * f_man)},
    "ATM C":       {"title": "Carbon in Atmosphere",        "col_idx": 4,  "unit": "(kg)",       "y2_label": r'$\mathrm{P_{atm}}$ (bar)', "y2_scale_factor": lambda y_val: y_val * g / (Aman * 1e5), "apply_conversion": True}, # ATM C に変換を適用
    "OC C":        {"title": "Carbon in Ocean",             "col_idx": 5,  "unit": "(kg)"},
    "P C":         {"title": "Carbon in Ocean Crust",       "col_idx": 6,  "unit": "(kg)"},
    "FD":          {"title": "Mid-Ocean Ridge Flux",        "col_idx": 10, "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    "FM":          {"title": "Arc Flux",                    "col_idx": 9,  "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    "Fplume":      {"title": "Plume Flux",                  "col_idx": 11, "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    "Ingas_UM":    {"title": "Ingassing Flux to Upper Mantle", "col_idx": 7, "unit": "(kg/yr)",   "ylim": (1e9, 1e15)},
    "Ingas_LM":    {"title": "Ingassing Flux to Lower Mantle", "col_idx": 8, "unit": "(kg/yr)",   "ylim": (1e9, 1e15)},
    "Inter":       {"title": "Lower to Upper Flux",         "col_idx": 12, "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    "Weathering_c": {"title": "Continental Weathering",     "col_idx": 14, "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    "Weathering_s": {"title": "Seafloor Weathering",        "col_idx": 15, "unit": "(kg/yr)",    "ylim": (1e9, 1e15)},
    # Melt_lowerの定義を修正: standardファイルではデータが存在しない、あるいはプロットしないため、col_idx_stdは不要
    "Melt_lower":  {"title": "Trapped Melt to Lower Mantle", "col_idx_100": 20, "col_idx_other": 18, "unit": "(kg/yr)", "ylim": (1e9, 1e15)}
}

# LaTeX表記マッピング (軸ラベル用)
xlabel_mapping = {
    "UM C": r'$\mathrm{M_{um}}$',
    "LM C": r'$\mathrm{M_{lm}}$',
    "M C": r'$\mathrm{M_{m}}$',
    "ATM C": r'$\mathrm{M_{atm}}$',
    "OC C": r'$\mathrm{M_{oc}}$',
    "P C": r'$\mathrm{M_{plate}}$',
    "FD": r'$\mathrm{F_{mor}}$',
    "FM": r'$\mathrm{F_{arc}}$',
    "Fplume": r'$\mathrm{F_{plume}}$',
    "Ingas_UM": r'$\mathrm{f_{slab}F_{ingas}}$',
    "Ingas_LM": r'$\mathrm{(1-f_{slab})F_{ingas}}$',
    "Inter": r'$\mathrm{F_{inter}}$',
    "Weathering_c": r'$\mathrm{F_{sil}}$',
    "Weathering_s": r'$\mathrm{F_{sw}}$',
    "Melt_lower": r'$\mathrm{f_{trap}F_{trap}}$'
}

# ファイルパスとf_trapの対応
DATA_FILES = {
    0: {"path": "carbon_median_standard.dat", "label": r'$\mathrm{f_{\mathrm{trap}}}=0$'},
    0.25: {"path": "carbon_median_25%.dat", "label": r'$\mathrm{f_{\mathrm{trap}}}=0.25$'},
    0.5: {"path": "carbon_median_50%.dat", "label": r'$\mathrm{f_{\mathrm{trap}}}=0.5$'},
    0.75: {"path": "carbon_median_75%.dat", "label": r'$\mathrm{f_{\mathrm{trap}}}=0.75$'},
    1.0: {"path": "carbon_median_100%.dat", "label": r'$\mathrm{f_{\mathrm{trap}}}=1$'},
}

# --- データ読み込み関数 ---
def load_median_data(file_path):
    """ 指定されたパスから全てのmedianデータを読み込み、辞書で返す """
    data_dict = {}
    if not os.path.exists(file_path):
        print(f"Warning: Data file not found: {file_path}. Skipping this dataset.")
        return None, None

    data = np.loadtxt(file_path)
    time = data[:, 0]

    for key, config in DATA_CONFIG.items():
        col_idx = None
        if key == "Melt_lower" or key == "Fres": # Melt_lower and Fres have special column indexing
            if "standard" in file_path:
                # Melt_lowerとFresはf_trap=0の場合はプロットしないため、ここではデータを読み込まない
                # そのため、col_idx_stdのようなキーは不要、またはここでNoneを設定
                col_idx = None 
            elif "100%" in file_path: # f_trap = 1.0 の場合
                col_idx = config.get("col_idx_100")
            else: # f_trap = 0.25, 0.5, 0.75 の場合
                col_idx = config.get("col_idx_other")
        else: # Standard column indexing for other variables
            col_idx = config.get("col_idx")

        if col_idx is not None and col_idx < data.shape[1]:
            data_dict[key] = data[:, col_idx]
        elif col_idx is not None:
            print(f"Warning: Column index {col_idx} for '{key}' out of bounds in {file_path}. Skipping.")
        # col_idx が None の場合は、その変数のデータは辞書に追加されない
            
    return time, data_dict

# 全てのデータセットをロード
all_times = {}
all_data_sets = {}
for f_trap_val, file_info in DATA_FILES.items():
    time, data_dict = load_median_data(file_info["path"])
    if time is not None:
        all_times[f_trap_val] = time
        all_data_sets[f_trap_val] = data_dict

# --- 複数プロットを1つのPDFにまとめる関数 ---
def generate_combined_plots(plot_keys, pdf_filename, plot_title_prefix=""):
    num_plots = len(plot_keys)
    ncols = 2 # 1行あたりの列数
    nrows = (num_plots + 1) // ncols # 必要な行数を計算 (切り上げ)

    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 22, nrows * 14)) # figsizeを調整
    axs = axs.flatten() # 1次元配列に変換

    plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.06, wspace=0.5, hspace=0.4)

    alphabet_iterator = iter(string.ascii_lowercase)

    for idx, key in enumerate(plot_keys):
        if idx >= len(axs):
            break

        ax1 = axs[idx]
        config = DATA_CONFIG.get(key)
        if not config:
            print(f"Error: Configuration for key '{key}' not found. Skipping plot.")
            ax1.text(0.5, 0.5, f'No config for {key}', ha='center', va='center', fontsize=25, color='gray')
            ax1.set_xticks([])
            ax1.set_yticks([])
            continue

        title = config["title"]
        unit = config["unit"]
        
        # Plot data for each f_trap value
        for f_trap_val in sorted(DATA_FILES.keys()):
            # Melt_lowerとFresについてはf_trap=0のデータをスキップ
            if (key == "Melt_lower" or key == "Fres") and f_trap_val == 0:
                continue # このf_trap値ではプロットしない

            if f_trap_val in all_data_sets and key in all_data_sets[f_trap_val]:
                time = all_times[f_trap_val]
                data = all_data_sets[f_trap_val][key]
                label = DATA_FILES[f_trap_val]["label"]
                
                # ATM C のデータに変換を適用
                if key == "ATM C" and config.get("apply_conversion"):
                    data = data * Aman * 1e5 / g

                ax1.plot(time, data, label=label, linewidth=2)
            else:
                pass # データがない場合はプロットしない

        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.tick_params(axis='x', which='major', labelsize=55, pad=15) # Adjust pad as needed
        ax1.tick_params(axis='y', which='major', labelsize=55)
        ax1.set_xlabel("Time (Gyr)", fontsize=65)
        ax1.set_ylabel(f"{xlabel_mapping.get(key, key)} {unit}", fontsize=65)
        ax1.set_title(title, fontsize=55)

        # y軸の範囲設定（'ylim'が設定されている場合のみ）
        if "ylim" in config:
            ylim_min, ylim_max = config["ylim"]
            ax1.set_ylim(ylim_min, ylim_max)
        
        ax1.legend(fontsize=30, loc="best")

        # 2nd y-axis (twinx) for specific variables
        if "y2_label" in config and "y2_scale_factor" in config:
            ax2 = ax1.twinx()
            
            # ax1のy軸範囲をax2のスケールファクタで変換してax2のy軸範囲を設定
            y_min, y_max = ax1.get_ylim()
            ax2_y_min = config["y2_scale_factor"](y_min)
            ax2_y_max = config["y2_scale_factor"](y_max)
            ax2.set_ylim(ax2_y_min, ax2_y_max)
                
            ax2.set_ylabel(config["y2_label"], fontsize=65)
            ax2.set_yscale("log")
            ax2.tick_params(axis='both', which='major', labelsize=55)
            
            # ax1とax2の凡例を統合
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            # 重複するラベルを削除し、ユニークなものだけを統合
            combined_labels = dict(zip(labels1, lines1))
            combined_labels.update(dict(zip(labels2, lines2)))
            ax1.legend(combined_labels.values(), combined_labels.keys(), fontsize=40, loc="best")

        else:
            # ax2がない場合はax1の凡例のみ表示
            ax1.legend(fontsize=40, loc="best")

        # アルファベットの追加
        try:
            label_char = next(alphabet_iterator)
            ax1.text(-0.25, 1.0, f'{label_char}', transform=ax1.transAxes,
                    fontsize=70, va='top', ha='left', weight='bold') 
        except StopIteration:
            pass

    # 残りの空白subplotを削除
    for j in range(num_plots, len(axs)):
        fig.delaxes(axs[j])

    # 1つのPDFファイルとして保存
    fig.savefig(pdf_filename, format='pdf', transparent=True, facecolor='white')
    plt.close(fig)
    print(f"✅ {plot_title_prefix} plots combined and saved to {pdf_filename}")

# --- メイン処理 ---
if __name__ == "__main__":
    # 炭素リザーバー関連のキーリスト (P C, ATM C の順序を調整)
    carbon_reservoir_keys = ["UM C", "LM C", "P C", "ATM C", "OC C"] 

    # 炭素フラックス関連のキーリスト
    carbon_flux_keys = ["FD", "FM", "Fplume", "Inter", "Ingas_UM", "Ingas_LM", "Weathering_c", "Weathering_s", "Melt_lower"]

    # ディレクトリの存在チェック（念のため）
    for f_trap_val, file_info in DATA_FILES.items():
        if not os.path.exists(file_info["path"]):
            print(f"Warning: Data file '{file_info['path']}' not found.")

    # 炭素リザーバーのプロットを生成
    print("\nGenerating Carbon Reservoir plots...")
    generate_combined_plots(carbon_reservoir_keys, "reservoir.pdf", "Carbon Reservoir")

    # 炭素フラックスのプロットを生成
    print("\nGenerating Carbon Flux plots...")
    generate_combined_plots(carbon_flux_keys, "flux.pdf", "Carbon Flux")

    print("\nAll plotting tasks completed.")