from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

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

def plot_performance_scatter(y_true, y_predict, eq_bounds=(30, 90), color='k', markersize=50):
    from sklearn.metrics import r2_score, mean_absolute_error
    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    ax.scatter(y_true, y_predict, color=color, edgecolors='#333333', s=markersize)
    ax.axline([eq_bounds[0], eq_bounds[0]], [eq_bounds[1], eq_bounds[1]], color='grey', ls='--')
    r, p = pearsonr(y_true, y_predict)
    ax.annotate(f'$r$ = {np.round(r, 3)}, P-val={"{0:.2e}".format(p)}', xy=[0.1, 0.9], xycoords='axes fraction')
    ax.annotate(f'$R^2$ = {np.round(r2_score(y_true, y_predict), 3)}', xy=[0.1, 0.85], xycoords='axes fraction')
    ax.annotate(f'MAE = {np.round(mean_absolute_error(y_true, y_predict), 2)}', xy=[0.1, 0.8], xycoords='axes fraction')
    ax.set_xlabel('True age, years')
    ax.set_ylabel('Predicted age, years')


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
                                                        'Predicted std':lambda x: np.sqrt(np.sum(np.square(x)))})
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
    fig.show()


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

