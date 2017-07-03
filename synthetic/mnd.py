from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import seaborn.apionly as sns
from scipy import stats


def fitMVN():
    # set plotting style
    sns.set_style("darkgrid")

    # load PM dataset
    pm_daily = np.loadtxt('../data/epa_hourly/alabamatest.csv', delimiter=',')
    #Nyears = int(np.shape(pm_daily)[0] / 365*24)
    Nsites = np.shape(pm_daily)[1]
    # Months = ['May', 'June', 'July', 'August', 'September', 'October', 'November', 'December', 'January', 'February',
    #           'March', 'April']


    # calculate standard deviation in daily flows each month and squared Mahalanobis distances
    #StdMonthly = calc_monthly_std(pm_daily, Nyears, Nsites)
    #D2 = calcD2(Nyears, Nsites, np.log(pm_daily))
    #D2 = calcD2(Nsites, np.log(pm_daily))

    # calculate theoretical quantiles for a chi^2 distribution with dof = Nsites, and for the standard normal distribution
    m = np.array(range(1, Nyears + 1))
    p = (m - 0.5) / Nyears
    chi2 = stats.chi2.ppf(p, Nsites)
    norm = stats.norm.ppf(p, 0, 1)

    # initialize matrices to store correlation coefficients and significance levels for marginal normal distributions and chi^2 distributions
    normCorr = np.zeros([Nsites, 12])
    norm_sigLevel = np.zeros([Nsites, 12])
    chi2Corr = np.zeros([12])
    chi2_sigLevel = np.zeros([12])

    for i in range(len(Months)):
        # plot histograms of standard deviation of daily flows each month, and of their logs
        plotHistograms(Nsites, StdMonthly[:, :, i], 'Standard Deviation of Daily ' + Months[i] + ' Flows',
                       Months[i] + 'Hist.png')
        plotHistograms(Nsites, np.log(StdMonthly[:, :, i]), 'log(Standard Deviation of Daily ' + Months[i] + ' Flows)', \
                       'Log' + Months[i] + 'Hist.png')

        # plot QQ plots of standard deviation of daily flows each month, and of their logs
        plotNormQQ(Nsites, StdMonthly[:, :, i], norm, 'Standard Deviation of Daily ' + Months[i] + ' Flows',
                   Months[i] + 'QQ.png')
        normCorr[:, i] = plotNormQQ(Nsites, np.log(StdMonthly[:, :, i]), norm,
                                    'log(Standard Deviation of Daily ' + Months[i] + ' Flows)',
                                    'Log' + Months[i] + 'QQ.png')

        # plot QQ plot of Chi Squared distribution of log of standard deviation in daily flows each month
        chi2Corr[i] = plotChi2QQ(Nsites, D2[:, i], chi2,
                                 'D$\mathregular{^2}\!$ of log(Standard Deviation of Daily ' + Months[i] + ' Flows)', \
                                 'Log' + Months[i] + 'Chi2QQ.png')

        # find significance levels
        chi2_sigLevel[i] = chi2_MC(Nsites, Nyears, chi2, chi2Corr[i])
        norm_sigLevel[:, i] = norm_MC(Nsites, Nyears, norm, normCorr[:, i])

    np.savetxt('Norm_sigLevels.txt', np.transpose(norm_sigLevel))
    np.savetxt('Norm_corr.txt', np.transpose(normCorr))
    np.savetxt('Chi2_sigLevels.txt', chi2_sigLevel)
    np.savetxt('Chi2_corr.txt', chi2Corr)

    return None


def calc_monthly_std(Qdaily, Nyears, Nsites):
    Nmonths = 12
    # first month = May (1st month of water year)
    DaysPerMonth = np.array([31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30])

    Qmonthly = np.zeros([Nsites, Nyears, Nmonths])
    StdMonthly = np.zeros([Nsites, Nyears, Nmonths])
    for year in range(Nyears):
        for month in range(Nmonths):
            start = year * 365 + np.sum(DaysPerMonth[0:month])

            for i in range(Nsites):
                # find total flow each month
                Qmonthly[i, year, month] = 86400 * np.sum(Qdaily[start:start + DaysPerMonth[month], i])

            # find standard deviation in daily flows each month
            for i in range(Nsites):
                for j in range(DaysPerMonth[month]):
                    StdMonthly[i, year, month] = StdMonthly[i, year, month] + \
                                                 (
                                                 86400 * Qdaily[start + j, i] - Qmonthly[i, year, month] / DaysPerMonth[
                                                     month]) ** 2

                StdMonthly[i, year, month] = np.sqrt((1 / (DaysPerMonth[month] - 1)) * StdMonthly[i, year, month])

    return StdMonthly


def plotHistograms(Nsites, data, xlabel, filename):
    fig = plt.figure()
    for i in range(Nsites):
        ax = fig.add_subplot(1, Nsites, i + 1)
        ax.hist(data[i, :], bins=10, color='navy', alpha=0.8)
        ax.set_title('Site ' + str(i + 1), fontsize=16)

    fig.text(0.1, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=14)
    fig.text(0.5, 0.04, xlabel, ha='center', fontsize=14)
    fig.subplots_adjust(bottom=0.15)
    fig.set_size_inches([22.525, 4.825])
    fig.savefig('Hists/' + filename)
    fig.clf()

    return None


def plotNormQQ(Nsites, data, norm, title, filename):
    corr = np.zeros([Nsites])
    fig = plt.figure()
    for i in range(Nsites):
        corr[i] = np.corrcoef(np.sort(data[i, :]), norm)[0, 1]
        z = (data[i, :] - np.mean(data[i, :])) / np.std(data[i, :])
        ax = fig.add_subplot(1, Nsites, i + 1)
        ax.scatter(norm, np.sort(z))
        ax.plot([-3, 3], [-3, 3], c='r')
        ax.set_title('Site ' + str(i + 1), fontsize=16)
        ax.set_xlim([-3, 3])
        ax.set_ylim([-3, 3])

    fig.text(0.1, 0.5, 'Sample Quantiles', va='center', rotation='vertical', fontsize=14)
    fig.text(0.5, 0.04, 'Theoretical Quantiles', ha='center', fontsize=14)
    fig.suptitle('Normal Q-Q Plot of ' + title, fontsize=16)
    fig.subplots_adjust(bottom=0.15, top=0.85)
    fig.set_size_inches([22.525, 4.825])
    fig.savefig('QQplots/' + filename)
    fig.clf()

    return corr


def calcD2(Nsites, data):
    D2 = np.zeros([Nyears, 12])
    X = np.zeros([Nyears, Nsites])
    Xprime = np.zeros([Nyears, Nsites])
    S = np.zeros(Nsites)
    for i in range(12):
        # fill data matrix, X, for ith month
        for j in range(Nsites):
            X[:, j] = data[j, :, i]

        # calculate covariance matrix, S, for ith month
        Xprime = X - (1 / Nyears) * np.dot(np.ones([Nyears, Nyears]), X)
        S = (1 / (Nyears - 1)) * np.dot(np.transpose(Xprime), Xprime)

        # calculate Mahalanobis distance, D2, for each year's ith month
        for j in range(Nyears):
            D2[j, i] = np.dot(np.dot((X[j, :] - np.mean(X, 0)), np.linalg.inv(S)),
                              (np.transpose(X[j, :] - np.mean(X, 0))))

    return D2


def plotChi2QQ(Nsites, data, chi2, title, filename):
    corr = np.corrcoef(np.sort(data), chi2)[0, 1]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(chi2, np.sort(data))
    ax.plot([0, 1.1 * np.max(chi2)], [0, 1.1 * np.max(chi2)], c='r')
    ax.set_xlabel('Theoretical Quantiles', fontsize=16)
    ax.set_xlim([0, 1.1 * np.max(chi2)])
    ax.set_ylabel('Sample Quantiles', fontsize=16)
    ax.set_ylim([0, 1.1 * np.max(data)])
    ax.tick_params(axis='both', labelsize=14)
    ax.set_title(r'$\chi^2$' + ' Q-Q Plot of ' + title, fontsize=16)
    fig.savefig('QQplots/' + filename)
    fig.clf()

    return corr


def chi2_MC(Nsites, Nyears, theoretical, dataCorr):
    corr = np.zeros(10000)
    for i in range(10000):  # 10,000 MC simulations
        simulated = stats.chi2.rvs(Nsites, size=Nyears)
        corr[i] = np.corrcoef(np.sort(simulated), theoretical)[0, 1]

    # find significance levels
    corr = np.sort(corr)
    for i in range(10000):
        if dataCorr > corr[i]:
            sigLevel = (i + 0.5) / 10000

    return sigLevel


def norm_MC(Nsites, Nyears, theoretical, dataCorr):
    sigLevel = np.zeros(Nsites)
    corr = np.zeros([10000])
    for i in range(10000):  # 10,000 MC simulations
        simulated = stats.norm.rvs(0, 1, size=Nyears)
        corr[i] = np.corrcoef(np.sort(simulated), theoretical)[0, 1]

    # find significance levels
    corr = np.sort(corr)
    for i in range(10000):
        for j in range(Nsites):
            if dataCorr[j] > corr[i]:
                sigLevel[j] = (i + 0.5) / 10000

    return sigLevel


fitMVN()