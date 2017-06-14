import pandas as pd
import matplotlib.pyplot as plt


SITE_NAMES = ['Bramblemet', 'Cambermet', 'Chimet', 'Sotonmet']
COL_NAMES = ["WSPD","WD","GST","ATMP","WTMP","BARO","DEPTH","AWVHT","MWVHT","WVHT","APD"]
COL_DESCRIPTIONS = ["Wind Speed","Wind Direction","Maximum Gust","Air Temperature","Water Temperature",
                    "Barometric Pressure","Water depth","Average Wave Height",
                    "Maximum Wave Height","Significant Wave Height","Average Wave Period"]
COL_UNITS = ["knots","degrees","knots","degrees C","degrees C","millibars","metres","metres",
             "metres","metres","seconds"]


def read_dataframe(index):
    # Read input file
    # Bramblemet
    df = pd.read_csv('../data/bramblemet/{}.csv'.format(SITE_NAMES[index]),
                     low_memory=False)

    # Parse Date Local
    df['Datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'],
                                    format='%d/%m/%Y %H:%M')
    df.set_index('Datetime', inplace=True)

    return df


def onresize(event):
    plt.tight_layout()


if __name__ == "__main__":

    siteId = 3;
    df = read_dataframe(siteId)

    for i in range(0,len(COL_NAMES)):
        if COL_NAMES[i] not in df.columns:
            continue

        if df[COL_NAMES[i]].isnull().sum() == len(df[COL_NAMES[i]]) :
            print('All the values of the column are NAN!')
            continue

        ts = pd.Series(df[COL_NAMES[i]])
        ts = ts.resample('5T').mean() # To fill missing dates/values

        fig = plt.figure()
        ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=2)
        #ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=2, wspace=0, hspace=0)
        ts.plot(ax=ax1)
        ax1.set_title(SITE_NAMES[siteId])
        ax1.set_ylabel("{}:{} [{}]".format(COL_NAMES[i],
                                       COL_DESCRIPTIONS[i],
                                       COL_UNITS[i]))

        ax2 = plt.subplot2grid((1, 3), (0, 2))
        ts.plot.hist(alpha=0.5, ax=ax2)
        # plt.subplots_adjust(left=0, bottom=0, right=0.01, top=0.01,
        #                     wspace=0, hspace=0)
        fig.tight_layout()
        fig.canvas.mpl_connect('resize_event', onresize)


    plt.show()
