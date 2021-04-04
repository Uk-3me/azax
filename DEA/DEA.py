from statsmodels.stats import multitest
from scipy import stats
from adjustText import adjust_text


def compare_clusters(orgX, comp, correction=True):
    # EMT group
    x_comp = pd.DataFrame(orgX.loc[:, X_label])

    # NonEMT group
    x_rest = pd.DataFrame(
        orgX.loc[:, [a for a in list(orgX) if a not in X_label]])

    # calculate mean
    comp_mean = np.mean(x_comp, axis=1)
    rest_mean = np.mean(x_rest, axis=1)

    # foldchange
    FCdown = rest_mean/comp_mean
    FCup = comp_mean/rest_mean

    # t-test p value
    P_V = []
    for a in tqdm(range(0, len(comp_mean))):
        stat, PV = stats.ttest_ind(
            x_comp.iloc[a, :], x_rest.iloc[a, :], equal_var=True)
        P_V.append(PV)

    # apply bonferroni correction
    if correction == True:
        rej, P_V = multitest.multipletests(P_V, method='b', alpha=0.05)[:2]

    # convert to df
    P_V = pd.DataFrame({'P_V': P_V}, index=FCup.index)
    FCup = pd.DataFrame({'FCup': FCup})
    FCdown = pd.DataFrame({'FCdown': FCdown})

    table_compup = pd.concat([FCup, P_V], axis=1)
    table_compdown = pd.concat([FCdown, P_V], axis=1)
    # sort
    table_compup = table_compup.sort_values('FCup', ascending=False)
    table_compdown = table_compdown.sort_values('FCdown', ascending=False)
    # if p-value<0.05, delete row
    table_filtup = table_compup[table_compup.P_V <= 0.05]
    table_filtdown = table_compdown[table_compdown.P_V <= 0.05]
    pd.set_option('display.max_rows', 200)
    return table_filtup


def volcano_plot(data, pvthresh, fcthresh):
    logFC = np.log2(data.FCup)
    logPV = -np.log10(list(data.P_V))
    genes = np.array([i for i in data.index])

    fig, ax = plt.subplots(figsize=[15, 10])
    sc = plt.scatter(logFC, logPV, s=10)
    ax.set_xlabel('log2(Fold Change)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title("Volcano plot of 5-aza treated EMT vs Non-EMT")

    FCthreshu = fcthresh
    PVthresh = pvthresh


#     plt.axvline(x=FCthreshu,ls='dotted')
#     plt.axvline(x=-FCthreshu, ls='dotted')
#     plt.axhline(y=PVthresh, ls='dotted')

    # only plot genes that has high FC and PV
    data2 = data.copy()
    data2.FCup = logFC
    data2.P_V = logPV
    textdfup = data2[(data2.FCup > fcthresh) & (data2.P_V > pvthresh)]
    textdfdown = data2[(data2.FCup < -fcthresh) & (data2.P_V > pvthresh)]
    updf = textdfup[0:20]
    concatdf = pd.concat([updf, textdfdown], axis=0)

    text = [plt.text(concatdf.FCup[i], concatdf.P_V[i],
                     concatdf.index[i], fontsize=10) for i in range(len(concatdf))]
    adjust_text(text)

    plt.show()
    return concatdf
