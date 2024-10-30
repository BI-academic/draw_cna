'''
    Date: 2024-03-27
    Author: duaghk
    Purpose: Draw chromosome scatterplot.
'''

# import library
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt



class PlotCNV:
    def __init__(self) -> None:
        # Set global value
        self.linestyle_dict = {
            'solid': 'solid',
            'dotted': 'dotted',
            'dashed': 'dashed',
            'dashdot': 'dashdot',
            'loosely dotted': (0, (1, 10)),
            'dotted': (0, (1, 1)),
            'densely dotted': (0, (1, 1)),
            'long dash with offset': (5, (10, 3)),
            'loosely dashed': (0, (5, 10)),
            'dashed': (0, (5, 5)),
            'densely dashed': (0, (5, 1)),
            'loosely dashdotted': (0, (3, 10, 1, 10)),
            'dashdotted': (0, (3, 5, 1, 5)),
            'densely dashdotted': (0, (3, 1, 1, 1)),
            'dashdotdotted': (0, (3, 5, 1, 5, 1, 5)),
            'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
            'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))
            }

    def draw_cna(
            self,
            data: pd.DataFrame, 
            x: str, 
            y: str, 
            chrom_list: list = None,
            y2: str = None,
            figsize: tuple = (24,3), 
            linestyle: str = "long dash with offset",
            palette_name: str = None,
            y2_color: str = None,
            ylim: tuple = (0, 400),
            purity: float = 1.0,
            title: str = None
        ):

        ref_chrom_list = [f"chr{x}" for x in range(1,23)] + ["chrX", "chrY"]
        ref_ratio_list = [5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 0.5]
        ref_ratio_dict = {k:v for k,v in zip(ref_chrom_list, ref_ratio_list)}
        if chrom_list is None:
            target_ratio_list = ref_ratio_list
            chrom_list = ref_chrom_list
        elif type(chrom_list) == list:
            target_ratio_list = [ref_ratio_dict[x] for x in chrom_list]
        elif type(chrom_list) == str:
            target_ratio_list = [ref_ratio_dict[chrom_list]]
            chrom_list = [chrom_list]

        # change datapoint to 2N
        data[y] *= 2
        if y2:
            data[y2] *= 2

        # set upper and lower threshold.
        ploidy = 2
        upper_limit = ploidy + (ploidy/2)*purity
        lower_limit = ploidy - (ploidy/2)*purity

        # set line style.
        linestyle_tuple = self.linestyle_dict[linestyle]
        # Get chromosome list.
        if palette_name: 
            color_pal = sns.color_palette(palette_name, len(chrom_list))
        else:
            # color_pal = ["black" if i % 2 == 0 else "gray" for i in range(len(chrom_list))]
            color_pal = ["#000000" if i % 2 == 0 else "#616161" for i in range(len(chrom_list))]
        if not y2_color:
            y2_color = "green"
            # y2_color = "#00c853"

        # Set figure with width ratios.        
        fig, axes = plt.subplots(
            ncols=len(chrom_list), 
            figsize=figsize,
            gridspec_kw={"width_ratios": target_ratio_list}
        )
        for i, chrom in enumerate(chrom_list):
            tmp_df = data[data['chrom'].isin([chrom])].copy().reset_index(drop=True)
            # Set break point.
            sns.scatterplot(
                data=tmp_df,
                x=x,
                y=y,
                alpha=0.2,
                s=5,
                ax=axes[i],
                color=color_pal[i]
            )
            if y2:
                # Split df for draw breaked.
                change_points = tmp_df.loc[tmp_df["change_point"]].index.tolist()
                change_points.append(tmp_df.index.max())
                prev_index = 0
                for idx in change_points:
                    seg_df = tmp_df.iloc[prev_index:idx]
                    # check segment value.
                    seg_value = seg_df[y2].unique().tolist()[0]
                    if seg_value <= lower_limit: 
                        # seg_color = "#2962ff"
                        seg_color = "blue"
                    elif seg_value >= upper_limit:
                        # seg_color = "#d50000"
                        seg_color = "red"
                    else:
                        seg_color = y2_color 
                    sns.lineplot(
                        data=tmp_df.iloc[prev_index:idx],
                        x=x,
                        y=y2,
                        alpha=0.7,
                        linewidth=2,
                        ax=axes[i],
                        color=seg_color
                    )
                    prev_index = idx + 1
            axes[i].set_xticks([tmp_df[x].median()], [chrom], rotation=90)
            axes[i].set_xlabel("")
            axes[i].set_xlim(tmp_df[x].min(), tmp_df[x].max())
            if i == 0:
                axes[i].set_ylabel("Ploidy")
            else:
                axes[i].set_ylabel("")
                axes[i].set_yticks([])
            axes[i].set_ylim(*ylim)
            axes[i].spines['right'].set(alpha=0.3, linestyle=linestyle_tuple)
            axes[i].spines['left'].set(alpha=0.3, linestyle=linestyle_tuple)
            # axes[i].spines['left'].set_visible(False)
        plt.subplots_adjust(wspace=0, hspace=0)
        if title:
            plt.suptitle(title)
        # plt.show()
        pass

    # def draw_target_chromosome(
    #         data: pd.DataFrame, 
    #         chrom_list: list,
    #         x: str, 
    #         y: str, 
    #         y2: str = None,
    #         figsize: tuple = (24,3), 
    #         linestyle: str = "long dash with offset",
    #         palette_name: str = None,
    #         y2_color: str = None,
    #         ylim: tuple = (0, 400),
    #         title: str = None
    #     ):
    #     # Set default value
    #     linestyle_dict = {
    #         'solid': 'solid',
    #         'dotted': 'dotted',
    #         'dashed': 'dashed',
    #         'dashdot': 'dashdot',
    #         'loosely dotted': (0, (1, 10)),
    #         'dotted': (0, (1, 1)),
    #         'densely dotted': (0, (1, 1)),
    #         'long dash with offset': (5, (10, 3)),
    #         'loosely dashed': (0, (5, 10)),
    #         'dashed': (0, (5, 5)),
    #         'densely dashed': (0, (5, 1)),
    #         'loosely dashdotted': (0, (3, 10, 1, 10)),
    #         'dashdotted': (0, (3, 5, 1, 5)),
    #         'densely dashdotted': (0, (3, 1, 1, 1)),
    #         'dashdotdotted': (0, (3, 5, 1, 5, 1, 5)),
    #         'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
    #         'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))
    #     }
    #     ref_chrom_list = [f"chr{x}" for x in range(1,23)] + ["chrX", "chrY"]
    #     ref_ratio_list = [5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 0.5]
    #     ref_ratio_dict = {k:v for k,v in zip(ref_chrom_list, ref_ratio_list)}
    #     target_ratio_list = [ref_ratio_dict[x] for x in chrom_list]

    #     # set line style.
    #     linestyle_tuple = linestyle_dict[linestyle]
    #     # Get chromosome list.
    #     if palette_name: 
    #         color_pal = sns.color_palette(palette_name, len(chrom_list))
    #     else:
    #         color_pal = ["black" if i % 2 == 0 else "gray" for i in range(len(chrom_list))]
    #     if not y2_color:
    #         y2_color = "green"
    #     # Calculate chromosome ratio for set plot width
    #     # chrom_len = []
    #     # for chrom in chrom_list:
    #     #     chrom_len.append(len(data[data["chrom"].isin([chrom])]))
    #     # chrom_len_2 = [round(x/min(chrom_len)) for x in chrom_len]
    #     # chromsome ratio is fixed. 

    #     # Set figure with width ratios.        
    #     fig, axes = plt.subplots(
    #         ncols=len(chrom_list), 
    #         figsize=figsize,
    #         gridspec_kw={"width_ratios": target_ratio_list}
    #     )
    #     for i, chrom in enumerate(chrom_list):
    #         tmp_df = data[data['chrom'].isin([chrom])]
    #         sns.scatterplot(
    #             data=tmp_df,
    #             x=x,
    #             y=y,
    #             alpha=0.1,
    #             s=5,
    #             ax=axes[i],
    #             color=color_pal[i]
    #         )
    #         if y2:
    #             # sns.scatterplot(
    #             sns.lineplot(
    #                 data=tmp_df,
    #                 x=x,
    #                 y=y2,
    #                 alpha=0.5,
    #                 linewidth=2,
    #                 # s=5,
    #                 ax=axes[i],
    #                 color=y2_color
    #             )
    #         axes[i].set_xticks([tmp_df[x].median()], [chrom], rotation=90)
    #         axes[i].set_xlabel("")
    #         axes[i].set_xlim(tmp_df[x].min(), tmp_df[x].max())
    #         if i != 0:
    #             axes[i].set_ylabel("")
    #             axes[i].set_yticks([])
    #         axes[i].set_ylim(*ylim)
    #         axes[i].spines['right'].set(alpha=0.3, linestyle=linestyle_tuple)
    #         axes[i].spines['left'].set(alpha=0.3, linestyle=linestyle_tuple)
    #         # axes[i].spines['left'].set_visible(False)
    #     plt.subplots_adjust(wspace=0, hspace=0)
    #     if title:
    #         plt.suptitle(title)
    #     # plt.show()
    #     pass





