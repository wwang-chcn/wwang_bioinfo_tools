from seaborn.categorical import _BoxPlotter


class box_plotter(_BoxPlotter):
    def add_statistic_annotate(self,
                               ax,
                               x1,
                               x2,
                               y,
                               h,
                               value,
                               text_type,
                               **kws):
        default_lw = mpl.rcParams["patch.linewidth"]
        lw = kws.get("linewidth", kws.get("lw", default_lw))
        width = x2 - x1
        x1, x2 = x1 + .1 * width, x2 - .1 * width
        ax.plot([x1, x1, x2, x2], [y + .5 * h, y + h, y + h, y + .5 * h], lw=1.5, c='k')
        if text_type == 'raw':
            ax.text((x1 + x2) * .5,
                    y + h * 1.5,
                    f'P = {value:.2e}',
                    ha='center',
                    va='bottom',
                    color='k')
        elif value >= 0.05:
            ax.text((x1 + x2) * .5,
                    y + h * 1.5,
                    'n.s.',
                    ha='center',
                    va='bottom',
                    color='k')
        elif value >= 0.01:
            ax.text((x1 + x2) * .5,
                    y + h * 1.5,
                    '*',
                    ha='center',
                    va='bottom',
                    color='k')
        elif value >= 0.001:
            ax.text((x1 + x2) * .5,
                    y_ + h * 1.5,
                    '**',
                    ha='center',
                    va='bottom',
                    color='k')
        else:
            ax.text((x1 + x2) * .5,
                    y + h * 1.5,
                    '***',
                    ha='center',
                    va='bottom',
                    color='k')
    def draw_statistic_annotate(
            self,
            ax,
            **kwargs):
        from operator import methodcaller
        from scipy import stats
        if kwargs.get('statistic_test') == 'mannwhitneyu':
            statistic_kwargs = {'use_continuity': kwargs.get('use_continuity', True), 
                                'alternative': kwargs.get('alternative', None)}
        elif kwargs.get('statistic_test') == 'wilcoxon':
            statistic_kwargs = {
                'zero_method': kwargs.get('zero_method', 'wilcox'),
                'correction': kwargs.get('correction', False),
                'alternative': kwargs.get('alternative', 'two-sided')
            }
        elif kwargs.get('statistic_test') == 'ttest_ind':
            statistic_kwargs = {'axis': kwargs.get('axis', 0), 'equal_var': kwargs.get('equal_var', True), 'nan_policy': kwargs.get('nan_policy', 'propagate')}
        elif kwargs.get('statistic_test') == 'ttest_rel':
            statistic_kwargs = {'axis': kwargs.get('axis', 0), 'nan_policy': kwargs.get('nan_policy', 'propagate')}
        else:
            raise ValueError(
                'statistic_test must be either mannwhitneyu, wilcoxon, ttest_ind or ttest_rel.'
            )
        ylim = ax.get_ylim()
        h = .05 * (ylim[1] - ylim[0])
        if self.plot_hues is None:
            for i in range(len(self.plot_data) - 1):
                pvalue = methodcaller(kwargs.get('statistic_test'),
                                      self.plot_data[i], self.plot_data[i + 1],
                                      **statistic_kwargs)(stats).pvalue
                y = max(mpl.cbook.boxplot_stats(self.plot_data[i]),
                        mpl.cbook.boxplot_stats(self.plot_data[i + 1]))
                self.add_statistic_annotate(
                    ax,
                    i,
                    i + 1,
                    y,
                    h,
                    pvalue,
                    text_type=kwargs.get('text_type','raw'))
        else:
            for i in range(len(self.plot_data)):
                offsets = self.hue_offsets
                for j in range(len(self.hue_names) - 1):
                    pvalue = methodcaller(
                        kwargs.get('statistic_test'), self.plot_data[i][
                            self.plot_hues[i] == self.hue_names[j]],
                        self.plot_data[i][self.plot_hues[i] == self.hue_names[
                            j + 1]], **statistic_kwargs)(stats).pvalue
                    y = max(
                        mpl.cbook.boxplot_stats(
                            self.plot_data[i][self.plot_hues[i] ==
                                              self.hue_names[j]])[0]['whishi'],
                        mpl.cbook.boxplot_stats(self.plot_data[i][
                            self.plot_hues[i] == self.hue_names[j + 1]])[0]
                        ['whishi']) + .25 * h
                    self.add_statistic_annotate(
                        ax,
                        i + offsets[j],
                        i + offsets[j + 1],
                        y,
                        h,
                        pvalue,
                        text_type=kwargs.get('text_type','raw'))
        ax.set_ylim(ylim[0], ylim[1] + 2.5 * h)


def boxplot(x=None,
            y=None,
            hue=None,
            data=None,
            order=None,
            hue_order=None,
            orient=None,
            color=None,
            palette=None,
            saturation=0.75,
            width=0.8,
            dodge=True,
            fliersize=5,
            linewidth=None,
            whis=1.5,
            ax=None,
            statistic_args = None,
            **kwargs):
    plotter = box_plotter(
        x,
        y,
        hue,
        data,
        order,
        hue_order,
        orient,
        color,
        palette,
        saturation,
        width,
        dodge,
        fliersize,
        linewidth)
    if ax is None:
        ax = plt.gca()
    kwargs.update(dict(whis=whis))
    plotter.draw_boxplot(ax, kwargs)
    plotter.annotate_axes(ax)
    if statistic_args:
        kwargs = statistic_args
        plotter.draw_statistic_annotate(ax,
                                    **kwargs)
    if plotter.orient == "h":
        ax.invert_yaxis()