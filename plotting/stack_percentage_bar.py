def stack_percentage_bar(results, horizontal=False):
    """\
    results : pandas.DataFrame
        results.index is the name of each sample.
        results.columns is the label of each category
    horizontal: bool
        If the figure should be horizontal.
    """
    
    # the following modules should been loaded before
    #import numpy as np
    #import matplotlib as mpl
    #mpl.rcParams['pdf.fonttype'] = 42
    #mpl.rcParams['ps.fonttype'] = 42
    #mpl.rcParams['font.sans-serif'] = 'Helvetica'
    #import matplotlib.pyplot as plt
    #import seaborn as sns
    
    percentage = results.div(results.sum(axis=1), axis=0) * 100
    percentage_cum = percentage.cumsum(axis=1)
    percentage_cum.loc[:, percentage_cum.shape[1] - 1] = 100.0
    category_colors = plt.get_cmap('RdYlGn')(np.linspace(
        0.15, 0.85, results.shape[1]))
    
    if horizontal:
        fig, ax = plt.subplots(figsize=(9.2, 5))
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
        
        for i, (category_name,
                color) in enumerate(zip(results.columns, category_colors)):
            widths = percentage.iloc[:, i]
            starts = percentage_cum.iloc[:, i] - widths
            ax.barh(results.index,
                    widths,
                    left=starts,
                    height=0.5,
                    label=category_name,
                    color=color)
            xcenters = starts + widths / 2
            
            r, g, b, _ = color
            text_color = 'white' if r * g * b < 0.33 else ('grey' if r * g * b < 0.66 else 'black')
            for y, (x, c) in enumerate(zip(xcenters, results.iloc[:, i])):
                ax.text(x,
                        y,
                        f'{c:,}',
                        ha='center',
                        va='center',
                        color=text_color)
        ax.set_xticks(np.arange(0, 120, 20))
        ax.set_xticklabels(np.arange(0, 120, 20))
        ax.set_xlabel('Percentage')
        ax.set_yticks(np.arange(percentage.shape[1] - 1))
        ax.set_yticklabels(percentage.index)
        ax.set_ylabel('Sample')
        ax.legend(ncol=results.shape[1],
                  bbox_to_anchor=(0, 1),
                  loc='lower left',
                  fontsize='small')
    else:
        fig, ax = plt.subplots(figsize=(8.4, 6.8))
        ax.set_ylim(0, 100)
        
        for i, (category_name,
                color) in enumerate(zip(results.columns, category_colors)):
            hights = percentage.iloc[:, i]
            starts = percentage_cum.iloc[:, i] - hights
            ax.bar(results.index,
                   hights,
                   bottom=starts,
                   width=0.5,
                   label=category_name,
                   color=color)
            xcenters = starts + hights / 2
            
            r, g, b, _ = color
            text_color = 'white' if r * g * b < 0.33 else ('grey' if r * g * b < 0.66 else 'black')
            for x, (c, y) in enumerate(zip(results.iloc[:, i], xcenters)):
                ax.text(x,
                        y,
                        f'{c:,}',
                        ha='center',
                        va='center',
                        color=text_color)
        ax.set_xticks(np.arange(percentage.shape[1] - 1))
        ax.set_xticklabels(percentage.index)
        ax.set_xlabel('Sample')
        ax.set_yticks(np.arange(0, 120, 20))
        ax.set_yticklabels(np.arange(0, 120, 20))
        ax.set_ylabel('Percentage')
        ax.legend(bbox_to_anchor=(1, 0), loc='lower left')
    return fig, ax
