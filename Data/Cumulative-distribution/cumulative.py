import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import string # アルファベット生成用
import os # ファイル存在チェック用

# 変数名リスト（表示するCDFの順番をここで指定）
variables = {
    "e_erupt": 4,
    "e_erupt_p": 5,
    "f_meta": 3,    
    "f_plume": 6,
    "slab_c": 2,
    "rate_mantle_c": 0,
    "part_c": 1
}

# LaTeX表記マッピング
xlabel_mapping = {
    "e_erupt": r'$\varepsilon_{\mathrm{erupt,mor}}$',
    "e_erupt_p": r'$\varepsilon_{\mathrm{erupt,p}}$',
    "f_meta": r'$\mathrm{f_{\mathrm{meta}}}$',
    "f_plume": r'$\mathrm{f_{\mathrm{plume}}}$',
    "rate_mantle_c": r'$\mathrm{f_{\mathrm{man}}}$',
    "part_c": r'$\mathrm{f_{\mathrm{part}}}$',
    "slab_c": r'$\mathrm{f_{\mathrm{slab}}}$'
}

# CDFデータファイルパスのテンプレート
cdf_data_paths = {
    "standard": "standard/{var_name}_cdf.dat",
    "melt-trap": "melt-trap/{var_name}_cdf.dat",
    "melt-trap-25": "melt-trap-25/{var_name}_cdf.dat",
    "melt-trap-50": "melt-trap-50/{var_name}_cdf.dat",
    "melt-trap-75": "melt-trap-75/{var_name}_cdf.dat",
}

def load_cdf_data(var_name, base_path_template):
    """CDFデータを読み込むヘルパー関数"""
    file_path = base_path_template.format(var_name=var_name)
    if os.path.exists(file_path):
        data = np.loadtxt(file_path)
        return data[:, 0], data[:, 1]
    else:
        print(f"Warning: CDF data file not found: {file_path}")
        return None, None

# --- 全てのCDFを1つのPDFにまとめる関数 ---
def generate_combined_cdfs(variables, xlabel_mapping, cdf_data_paths):
    num_plots = len(variables)
    ncols = 2 # 1行あたりの列数
    nrows = (num_plots + 1) // ncols # 必要な行数を計算 (切り上げ)

    # FigureとAxesの作成
    # figsizeを調整して、全体のバランスを良くする
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 16, nrows * 12)) 
    axs = axs.flatten() # 1次元配列に変換

    # サブプロット間の余白を調整
    plt.subplots_adjust(left=0.1, right=0.95, top=0.97, bottom=0.05, wspace=0.3, hspace=0.3)

    # アルファベット生成用のイテレータ
    alphabet_iterator = iter(string.ascii_lowercase)

    for idx, var_name in enumerate(variables):
        if idx >= len(axs): # 念のため、サブプロットの数を超えないようにチェック
            break

        ax = axs[idx]
        xlabel = xlabel_mapping.get(var_name, var_name)
        
        # 各バージョンのCDFデータを読み込み
        x_standard, y_standard = load_cdf_data(var_name, cdf_data_paths["standard"])
        x_mt, y_mt = load_cdf_data(var_name, cdf_data_paths["melt-trap"])
        x_mt25, y_mt25 = load_cdf_data(var_name, cdf_data_paths["melt-trap-25"])
        x_mt50, y_mt50 = load_cdf_data(var_name, cdf_data_paths["melt-trap-50"])
        x_mt75, y_mt75 = load_cdf_data(var_name, cdf_data_paths["melt-trap-75"])

        # データが存在する場合のみプロット
        if x_standard is not None:
            ax.plot(x_standard, y_standard, lw=5, label=f"{r'$\mathrm{f_{\mathrm{trap}}}=0$'}")
        if x_mt25 is not None:
            ax.plot(x_mt25, y_mt25, lw=5, label=f"{r'$\mathrm{f_{\mathrm{trap}}}=0.25$'}")
        if x_mt50 is not None:
            ax.plot(x_mt50, y_mt50, lw=5, label=f"{r'$\mathrm{f_{\mathrm{trap}}}=0.5$'}")
        if x_mt75 is not None:
            ax.plot(x_mt75, y_mt75, lw=5, label=f"{r'$\mathrm{f_{\mathrm{trap}}}=0.75$'}")
        if x_mt is not None:
            ax.plot(x_mt, y_mt, lw=5, label=f"{r'$\mathrm{f_{\mathrm{trap}}}=1$'}")

        ax.grid(True, linestyle='--', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=45)
        ax.set_xlabel(xlabel, fontsize=55)
        ax.set_ylabel('Cumulative Probability', fontsize=55)

        # x軸上限の設定 (変数ごとに異なる場合)
        x_max = 0.5 if var_name == "f_plume" else 1.0
        ax.set_xlim(0, x_max)
        ax.set_ylim(-0.1, 1.1)
        ax.legend(fontsize=35)

        # アルファベットの追加
        try:
            label_char = next(alphabet_iterator)
            ax.text(-0.24, 1.0, f'{label_char}', transform=ax.transAxes,
                    fontsize=50, va='top', ha='left', weight='bold') 
        except StopIteration:
            pass

    # 残りの空白subplotを削除
    for j in range(num_plots, len(axs)):
        fig.delaxes(axs[j])

    # 1つのPDFファイルとして保存
    output_pdf_file = 'cdf.pdf'
    fig.savefig(output_pdf_file, format='pdf', transparent=True, facecolor='white')
    plt.close(fig)
    print(f"✅ All CDFs combined and saved to {output_pdf_file}")

# --- メイン処理 ---
if __name__ == "__main__":
    # 必要なディレクトリの存在チェック
    required_dirs = ["standard", "melt-trap", "melt-trap-25", "melt-trap-50", "melt-trap-75"]
    for d in required_dirs:
        if not os.path.exists(d):
            print(f"Warning: Directory '{d}' not found. Please ensure data files are in correct subdirectories.")

    print("Generating combined CDF plots...")
    generate_combined_cdfs(variables, xlabel_mapping, cdf_data_paths)
    print("\nAll plotting tasks completed.")