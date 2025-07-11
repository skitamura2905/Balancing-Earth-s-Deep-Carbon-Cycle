import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib.backends.backend_pdf import PdfPages
import string # アルファベット生成用
import os # ファイル存在チェック用

# --- 定数と設定 ---
# ヒストグラム用のファイル範囲
original_hist_file_range = range(2)  # original データ
new_hist_file_range = range(1, 10001)  # 新しいデータ (parameter)

# 読み込み対象の変数と列インデックスの対応表
variables = {
    "e_erupt": 2,
    "f_meta": 1,
    "rate_mantle_c": 0
}

# ラベルのLaTeX表記マッピング (共通)
label_mapping = {
    "e_erupt": r'$\varepsilon_{\mathrm{erupt,mor}}$',
    "e_erupt_p": r'$\varepsilon_{\mathrm{erupt,p}}$',
    "f_meta": r'$\mathrm{f_{\mathrm{meta}}}$',
    "f_plume": r'$\mathrm{f_{\mathrm{plume}}}$',
    "rate_mantle_c": r'$\mathrm{f_{\mathrm{man}}}$',
    "part_c": r'$\mathrm{f_{\mathrm{part}}}$',
    "slab_c": r'$\mathrm{f_{\mathrm{slab}}}$'
}


# --- データ読み込み関数 ---
def load_data_for_variable(base_dir, file_range, column_index, is_original=False):
    """
    指定されたディレクトリ、ファイル範囲、列インデックスに基づいてデータを読み込む。
    is_original=True の場合、データは [:, column_index] で読み込まれる。
    """
    all_values = []
    for i in file_range:
        # ディレクトリパスを正しく結合
        file_path = os.path.join(base_dir, f'parameter_{i}.dat')
        try:
            if not os.path.exists(file_path):
                continue
            
            data = np.loadtxt(file_path)
            if is_original: # original データの場合、2D配列の可能性があるので[:, column_index]
                # dataが1次元配列で[:, column_index]がIndexErrorになるのを防ぐ
                if data.ndim == 1:
                    values = np.atleast_1d(data[column_index])
                else:
                    values = np.atleast_1d(data[:, column_index])
            else: # parameter データの場合、1D配列の可能性があるので[column_index]
                values = np.atleast_1d(data[column_index])
            
            values = values[~np.isnan(values)] # NaN を除去
            all_values.extend(values)
        except IndexError:
            # print(f"Warning: Column index {column_index} out of bounds for {file_path}. Skipping.")
            continue
        except Exception as e:
            # print(f"An error occurred while reading {file_path}: {e}. Skipping.")
            continue
    return all_values


# --- 複数ヒストグラムを1つのPDFにまとめる関数 (単一データセット用) ---
def generate_combined_histograms(variables, label_mapping, file_range):
    num_histograms = len(variables)
    ncols = 3 # 1行あたりの列数
    nrows = (num_histograms + ncols - 1) // ncols # 必要な行数を計算 (切り上げ)

    # FigureとAxesの作成
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 10, nrows * 9)) # figsizeを調整
    axs = axs.flatten() # 1次元配列に変換
 
    # サブプロット間の余白を調整
    plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.16, wspace=0.3, hspace=0.4)

    # アルファベット生成用のイテレータ
    alphabet_iterator = iter(string.ascii_lowercase)

    for idx, (var_name, column_index) in enumerate(variables.items()):
        if idx >= len(axs): # 念のため、サブプロットの数を超えないようにチェック
            break

        ax = axs[idx]
        # === ここからインデントを修正 ===
        all_values = load_data_for_variable('parameter', file_range, column_index) # parameterディレクトリから読み込み

        xlabel_name = label_mapping.get(var_name, var_name) # ループ内で定義

        if all_values:
            ax.hist(all_values, bins=30, color='blue', edgecolor='black', alpha=0.7)
            ax.tick_params(axis='both', which='major', labelsize=40)
            
            ax.set_xlabel(xlabel_name, fontsize=50)
            ax.set_ylabel('Frequency', fontsize=50)
            #ax.set_title(f'Histogram of {xlabel_name}', fontsize=32) # タイトルはコメントアウトのまま
        else:
            ax.text(0.5, 0.5, f'No data for {xlabel_name}', ha='center', va='center', fontsize=25, color='gray')
            ax.set_xlabel(xlabel_name, fontsize=50)
            ax.set_ylabel('Frequency', fontsize=50)
            #ax.set_title(f'Histogram of {xlabel_name}', fontsize=32)
            ax.set_xticks([])
            ax.set_yticks([])

        # アルファベットの追加
        try:
            label_char = next(alphabet_iterator)
            ax.text(-0.25, 1.0, f'{label_char}', transform=ax.transAxes,
                    fontsize=45, va='top', ha='left', weight='bold') # フォントサイズと太さを調整
        except StopIteration:
            pass
        # === ここまでインデントを修正 ===

    # 残りの空白subplotを削除
    for j in range(num_histograms, len(axs)):
        fig.delaxes(axs[j])

    # 1つのPDFファイルとして保存
    output_pdf_file = 'hist.pdf'
    fig.savefig(output_pdf_file, format='pdf', transparent=True, facecolor='white')
    plt.close(fig)
    print(f"✅ All histograms combined and saved to {output_pdf_file}")


# --- 比較ヒストグラムを1つのPDFにまとめる関数 ---
def generate_combined_comparison_histograms(variables, label_mapping, 
                                            original_file_range, new_file_range):
    num_histograms = len(variables)
    ncols = 3 
    nrows = (num_histograms + ncols - 1) // ncols

    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 10, nrows * 9)) 
    axs = axs.flatten()

    plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.16, wspace=0.3, hspace=0.4)

    alphabet_iterator = iter(string.ascii_lowercase)

    for idx, (var_name, column_index) in enumerate(variables.items()):
        if idx >= len(axs):
            break

        ax = axs[idx]
        
        # データの読み込み
        # load_data_for_variableのbase_dir引数を正しく指定
        original_values = load_data_for_variable('original', original_file_range, column_index, is_original=True)
        new_values = load_data_for_variable('parameter', new_file_range, column_index)

        xlabel_name = label_mapping.get(var_name, var_name)

        if original_values or new_values:
            all_values_for_bins = []
            if original_values: all_values_for_bins.extend(original_values)
            if new_values: all_values_for_bins.extend(new_values)
            
            bins = np.histogram_bin_edges(all_values_for_bins, bins=30)

            # original データのヒストグラム
            if original_values:
                ax.hist(original_values, bins=bins, color='blue', edgecolor='black', alpha=0.5, label='Original Data')

            # 新しいデータのヒストグラム
            if new_values:
                ax.hist(new_values, bins=bins, color='red', edgecolor='black', alpha=0.5, label='Successful Data')

            ax.tick_params(axis='both', which='major', labelsize=40)
            ax.set_xlabel(xlabel_name, fontsize=50)
            ax.set_ylabel('Frequency', fontsize=50)
            ax.set_yscale("log") # Y軸を対数スケールに
            ax.legend(fontsize=30) # 凡例を表示
            #ax.grid(True) # グリッドを追加
            #ax.set_ylim(1,1e5) # Y軸の範囲を固定したい場合はコメント解除
        else:
            ax.text(0.5, 0.5, f'No data for {xlabel_name}', ha='center', va='center', fontsize=25, color='gray')
            ax.set_xlabel(xlabel_name, fontsize=50)
            ax.set_ylabel('Frequency', fontsize=50)
            #ax.set_title(f'Histogram of {xlabel_name}', fontsize=32) # タイトルは表示するか検討
            ax.set_xticks([])
            ax.set_yticks([])

        # アルファベットの追加
        try:
            label_char = next(alphabet_iterator)
            ax.text(-0.25, 1.0, f'{label_char}', transform=ax.transAxes,
                    fontsize=45, va='top', ha='left', weight='bold') 
        except StopIteration:
            pass

    # 残りの空白subplotを削除
    for j in range(num_histograms, len(axs)):
        fig.delaxes(axs[j])

    output_pdf_file = 'hist_compare.pdf'
    fig.savefig(output_pdf_file, format='pdf', transparent=True, facecolor='white')
    plt.close(fig)
    print(f"✅ All comparison histograms combined and saved to {output_pdf_file}")


# --- メイン処理 ---
if __name__ == "__main__":
    # ディレクトリの存在チェック
    if not os.path.exists('parameter'):
        print("Warning: 'parameter' directory not found. Please ensure your data files are in a 'parameter' subdirectory.")
    if not os.path.exists('original'):
        print("Warning: 'original' directory not found. Please ensure your data files are in an 'original' subdirectory.")
    
    # 単一データセットヒストグラムの生成
    print("Generating combined histograms (single dataset)...")
    generate_combined_histograms(variables, label_mapping, new_hist_file_range) # parameterデータを使用

    # 比較ヒストグラムの生成
    print("\nGenerating combined comparison histograms...")
    generate_combined_comparison_histograms(variables, label_mapping, 
                                            original_hist_file_range, new_hist_file_range)

    print("\nAll plotting tasks completed.")