from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
import warnings
from sklearn.decomposition import PCA
from .shift import calculate_shift

METHYLCLOCK_NAMES_MAPPER = {
    'Horvath':'Horvath (353 CpGs)', 
    'Hannum':'Hannum (71 CpGs)', 
    'Levine':'Levine (513 CpGs)', 
    'BNN':'BNN (353 CpGs)', 
    'skinHorvath': 'Horvath skin (391 CpGs)', 
    'PedBE': 'PedBE (84 CpGs)', 
    'Wu': 'Wu (111 CpGs)',
    'TL': 'TL (140 CpGs)', 
    'BLUP': 'BLUP (319607 CpGs)', 
    'EN': 'Zhang (514 CpGs)', 
    }

def plot_performance_scatter(y_true, y_predict, eq_bounds=(30, 90), color='k', markersize=50, units='years'):
    from sklearn.metrics import r2_score, mean_absolute_error
    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    ax.scatter(y_true, y_predict, color=color, edgecolors='#333333', s=markersize)
    ax.axline([eq_bounds[0], eq_bounds[0]], [eq_bounds[1], eq_bounds[1]], color='grey', ls='--')
    r, p = pearsonr(y_true, y_predict)
    ax.annotate(f'$r$ = {np.round(r, 3)}, P-val={"{0:.2e}".format(p)}', xy=[0.1, 0.9], xycoords='axes fraction')
    ax.annotate(f'$R^2$ = {np.round(r2_score(y_true, y_predict), 3)}', xy=[0.1, 0.85], xycoords='axes fraction')
    ax.annotate(f'MAE = {np.round(mean_absolute_error(y_true, y_predict), 2)}, {units}', xy=[0.1, 0.8], xycoords='axes fraction')
    ax.set_xlabel(f'True age, {units}')
    ax.set_ylabel(f'Predicted age, {units}')


def plot_repr_uncertainty(y_day, y_predict, y_predict_std, 
                          days=None, 
                          nstd=2, 
                          dh=.05, 
                          barh=.05,
                          xlabel='Day of reprogramming',
                          ylabel='Predicted age, years',
                          ):
    print('Avg uncertainty std:', y_predict_std.mean())
    pred_df = pd.DataFrame({'Reprogramming day':y_day, 
                            'Predicted age':y_predict,
                            'Predicted std':y_predict_std,
                            })
    pred_df = pred_df.groupby('Reprogramming day').agg({'Predicted age':np.mean, 
                                                        'Predicted std':lambda x: np.sqrt(np.mean(np.square(x)))})
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.scatter(y_day, y_predict)
    ax.plot(pred_df.index, pred_df['Predicted age'], label='Average predicted age')
    ax.fill_between(pred_df.index, 
                 pred_df['Predicted age'] + nstd * pred_df['Predicted std'], 
                 pred_df['Predicted age'] - nstd * pred_df['Predicted std'], 
                 alpha=0.1, color='green', label='Credible interval')
    ax.axhline(0, color='grey', ls='--')
    ax.set_xticks(pred_df.index)
    ax.set_xlim([min(pred_df.index) - 0.5, max(pred_df.index) + 0.5])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    
    if days is not None:
        if type(dh) is not list:
            dh = [dh] * len(days)
        for d, hh in zip(days, dh):
            from pymare import meta_regression
            cond = (y_day == d[0]) | (y_day == d[1])
            y_rep_pred_mean_sample = y_predict[cond]
            y_rep_pred_std_sample = y_predict_std[cond]
            y_rep_day_sample = y_day[cond]

            # print(y_rep_day_sample)
            # print(y_rep_pred_mean_sample)
            # print(y_rep_pred_std_sample)

            meta_df_ = meta_regression(y_rep_pred_mean_sample, 
                            y_rep_pred_std_sample**2, 
                            (y_rep_day_sample == d[0]).astype(int),
                            X_names=['rej'], 
                            add_intercept=True, 
                            method='REML').to_df().set_index('name')
            pval = meta_df_.loc['rej', 'p-value']
            print(f'P-value of rejuvenation effect between {d[0]} and {d[1]} days is {pval}')
            barplot_annotate_brackets(d[0], d[1], 
                                    pval, pred_df.index, pred_df['Predicted age'],
                                    dh=hh, barh=barh)


#taken from https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    # lx, ly = center.get_loc(num1), height.loc[num1]
    # rx, ry = center.get_loc(num2), height.loc[num2]
    lx, ly = num1, height.loc[num1]
    rx, ry = num2, height.loc[num2]

    # if yerr:
    #     ly += yerr[num1]
    #     ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = min(ly, ry) - dh

    barx = [lx, lx, rx, rx]
    bary = [y, y-barh, y-barh, y]
    mid = ((lx+rx)/2, y-2*barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='center')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)

def plot_covariate_panel(Xa_clock, Xr_clock, ya, yr,
                        cbar_color_train = 'Blues',
                        cbar_color_test = 'Reds',
                        color_train = '#86d2da',
                        color_test = '#f24734',
                        xa_thr = 0.05,
                        xr_thr = 0.05,
                        cbar_label_train = 'Age, years',
                        cbar_label_test = 'Reprogramming day',
                        legend_label_train = 'Aging',
                        legend_label_test = 'Reprogramming',
                        ks_hist_title = 'Histogram of Kolmogorov-Smirnov test\n BH-adjusted p-values for Horvath skin sites',
                        pca_title = '',
                        pca_legend_pos = 'lower right',
                        pca_xlim = None,
                        pca_ylim = None,
                        custom_bins = None,
                        ):

    if custom_bins is not None:
        hist_bins = custom_bins 
    else:
        hist_bins = 16
    warnings.filterwarnings('ignore')

    #conduct PCA on joined dataset
    pca = PCA(2)
    pca.fit(pd.concat([Xa_clock, Xr_clock]))
    Xpa = pca.transform(Xa_clock)
    Xpr = pca.transform(Xr_clock)
    pvar = np.round(pca.explained_variance_ratio_, 2)

    #plotting starts here
    fig = plt.figure(figsize=(18, 5), constrained_layout=True)
    gs = GridSpec(2, 7, figure=fig, width_ratios=[4, 4, 0.5, 3, 3, 4, 4], wspace=0.8, hspace=0.5) 

    axcbrs = gs[:, 2].subgridspec(2, 1)
    axcb1, axcb2 = axcbrs.subplots()

    minihists = gs[:, 3:5].subgridspec(2, 2, wspace=0.3, hspace=0.5)
    minihist_subs = minihists.subplots(sharex=True)
    ax0 = fig.add_subplot(gs[:, :2])
    ax1 = minihist_subs[0,0]
    ax2 = minihist_subs[0,1]
    ax3 = minihist_subs[1,0]
    ax4 = minihist_subs[1,1]
    ax5 = fig.add_subplot(gs[:, 5:])

    #cbars
    norm_age = mpl.colors.Normalize(vmin=min(ya), vmax=max(ya))
    norm_rep = mpl.colors.Normalize(vmin=min(yr), vmax=max(yr))
    cb_age = mpl.colorbar.ColorbarBase(axcb2, cmap=cbar_color_train, norm=norm_age)
    cb_rep = mpl.colorbar.ColorbarBase(axcb1, cmap=cbar_color_test, norm=norm_rep)
    cb_age.set_label(cbar_label_train, loc='bottom')
    cb_rep.set_label(cbar_label_test, loc='bottom')
    axcb1.yaxis.set_ticks_position('left')
    axcb2.yaxis.set_ticks_position('left')

    # ax0 - PCA
    ax0.set_title(pca_title)
    sns.kdeplot(Xpa[:, 0], Xpa[:, 1], ax=ax0, color=color_train, levels=5, thresh=xa_thr)
    sns.kdeplot(Xpr[:, 0], Xpr[:, 1], ax=ax0, color=color_test, levels=5, thresh=xr_thr)
    z1_plot = ax0.scatter(Xpa[:, 0], Xpa[:, 1], c=ya, s=50, label=legend_label_train, 
                        edgecolors='#333333', cmap=cbar_color_train, zorder=4)
    z2_plot = ax0.scatter(Xpr[:, 0], Xpr[:, 1], c=yr, marker='s', s=50, label=legend_label_test, 
                        edgecolors='#333333', cmap=cbar_color_test, zorder=4)
    ax0.set_xlabel(f'PCA1: {int(pvar[0]*100)}%')
    ax0.set_ylabel(f'PCA2: {int(pvar[1]*100)}%')
    ax0.legend(loc=pca_legend_pos)
    leg = ax0.get_legend()
    leg.legendHandles[0].set(edgecolor='#333333')
    leg.legendHandles[1].set(edgecolor='#333333')
    ax0.set_xlim(pca_xlim)
    ax0.set_ylim(pca_ylim)

    # ax5 - overall KS pvalues hist
    ks, ksp = calculate_shift(np.array(Xa_clock), np.array(Xr_clock), metric_type='ks')
    kspa = multipletests(ksp, method='hs')[1]
    kspa = pd.DataFrame(index=Xa_clock.columns, data={'pval':kspa})

    kspaval = -np.log10(kspa['pval'].values)
    kspaval = np.where(kspaval < 100, kspaval , 100)
    pct = (kspaval < 2).sum() / kspaval.shape[0] * 100
    print(f'{round(pct, 2)}% of sites are not rejected by KS test')

    ax5.set_title(ks_hist_title)
    ax5.hist(kspaval, density=False, alpha=0.7, color='grey', bins=hist_bins)
    ax5.axvline(-np.log10(0.01), color='k', ls='--', lw=2.0)
    ax5.set_xlabel('KS test, -log10(adj.P-value)')
    ax5.set_ylabel('Number of CpGs')

    #top sites from Skin clock based on correlation with age
    selected_sites = Xa_clock.corrwith(ya).abs().sort_values(ascending=False).index[:4]

    #examples of sites + PSIs
    fig.text(0.54, 0.04, 'CpG methylation beta balue', ha='center')
    for k, ax_ in enumerate([ax1, ax2, ax3, ax4]): 
        site_pval = kspa.loc[selected_sites[k], 'pval']
        ax_.hist(Xa_clock[selected_sites[k]], density=True, alpha=0.8, 
                color=color_train, 
                bins=np.linspace(0, 1, 25), 
                edgecolor='none')
        ax_.hist(Xr_clock[selected_sites[k]], density=True, alpha=0.7, 
                color=color_test, 
                bins=np.linspace(0, 1, 25), 
                edgecolor='none')
        ax_.set_title(f'{selected_sites[k]},\n P-val={"{0:.2e}".format(site_pval)}')
        ax_.set_xlim([-0.05, 1.05])