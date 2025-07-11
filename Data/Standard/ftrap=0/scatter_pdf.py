import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib.backends.backend_pdf import PdfPages
import string # アルファベット生成用

# ファイル範囲
file_range = range(1, 10001)

# 読み込み対象の変数と列インデックスの対応表
variables = {
    "e_erupt": 4,
    "e_erupt_p": 5,
    "f_meta": 3,
    "f_plume": 6,
    "slab_c": 2,
    "rate_mantle_c": 0,
    "part_c": 1
}

# ラベルのLaTeX表記マッピング
label_mapping = {
    "e_erupt": r'$\varepsilon_{\mathrm{erupt,mor}}$',
    "e_erupt_p": r'$\varepsilon_{\mathrm{erupt,p}}$',
    "f_meta": r'$\mathrm{f_{\mathrm{meta}}}$',
    "f_plume": r'$\mathrm{f_{\mathrm{plume}}}$',
    "rate_mantle_c": r'$\mathrm{f_{\mathrm{man}}}$',
    "part_c": r'$\mathrm{f_{\mathrm{part}}}$',
    "slab_c": r'$\mathrm{f_{\mathrm{slab}}}$'
}

# 変数の組み合わせ (7C2 通り)
variable_pairs = list(itertools.combinations(variables.keys(), 2))

# 出力用PDF
# rasterized_dpi: ラスタライズするときのDPIを設定（適宜調整してください）
# たとえば、600 は高品質ですがファイルサイズは大きくなります。
# 300 はバランスの取れた選択肢です。
with PdfPages("scatter.pdf", metadata={'CreationDate': None}) as pdf: # ファイル名を変更
    fig, axs = plt.subplots(nrows=7, ncols=3, figsize=(21, 45))  # 3列構成, 21個で7行必要
    axs = axs.flatten()

    # === subplots_adjust で余白を調整 ===
    plt.subplots_adjust(left=0.08, right=0.95, top=0.98, bottom=0.03, wspace=0.5, hspace=0.5)

    # アルファベット生成用のイテレータ
    alphabet_iterator = iter(string.ascii_lowercase)

    for idx, (var_x, var_y) in enumerate(variable_pairs):
        # 変数ペアの数がサブプロットの数を超えたらループを抜ける
        if idx >= len(axs):
            break

        all_x = []
        all_y = []

        index_x = variables[var_x]
        index_y = variables[var_y]

        # データ取得
        for i in file_range:
            file_path = f'parameter/parameter_{i}.dat'
            try:
                data = np.loadtxt(file_path)
                x_values = np.atleast_1d(data[index_x])
                y_values = np.atleast_1d(data[index_y])
                all_x.extend(x_values)
                all_y.extend(y_values)
            except FileNotFoundError:
                continue
            except IndexError:
                continue
            except Exception as e:
                continue

        # プロット
        ax = axs[idx]
        if all_x and all_y:
            # === ここに rasterized=True を追加 ===
            # 大量の点を扱う際に、これをTrueにすることでPDFのファイルサイズを削減できます。
            # また、dpiを指定することで、ラスタライズされた画像の解像度を調整できます。
            ax.scatter(all_x, all_y, alpha=0.5, color='blue', edgecolors='black', linewidth=0.2, rasterized=True)
            ax.set_xlabel(label_mapping.get(var_x, var_x), fontsize=50)
            ax.set_ylabel(label_mapping.get(var_y, var_y), fontsize=50)
            ax.tick_params(axis='both', which='major', labelsize=35)
            ax.grid(True)
        else:
            ax.text(0.5, 0.5, 'No data available for this pair', ha='center', va='center', fontsize=30, color='gray')
            ax.set_xlabel(label_mapping.get(var_x, var_x), fontsize=50)
            ax.set_ylabel(label_mapping.get(var_y, var_y), fontsize=50)
            ax.tick_params(axis='both', which='major', labelsize=35)
            ax.set_xticks([])
            ax.set_yticks([])

        # === ここにアルファベットの追加 ===
        try:
            label_char = next(alphabet_iterator)
            ax.text(-0.35, 1.0, f'{label_char}', transform=ax.transAxes,
                    fontsize=45, va='top', ha='left', weight='bold') # フォントサイズと太さを調整
        except StopIteration:
            # アルファベットが尽きた場合（26個以上サブプロットがある場合）
            pass


    # 残りの空白subplotを消す
    for j in range(len(variable_pairs), len(axs)):
        fig.delaxes(axs[j])

    # pdf.savefig() に dpi を指定して、ラスタライズされたコンテンツの解像度を設定
    pdf.savefig(fig, dpi=300) # ここでDPIを設定
    plt.close()

print("✅ scatter.pdf を生成しました。")