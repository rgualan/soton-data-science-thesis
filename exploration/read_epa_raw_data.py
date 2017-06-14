import pandas as pd
import matplotlib.pyplot as plt


def read_dataframe():
    # Read input file
    df = pd.read_csv('../data/epa_raw/Alabama/AQDM_907365160.txt', low_memory=False)

    # Parse Date Local
    df['datetime'] = pd.to_datetime(df['datetime'],
                                    format='%Y%m%dT%H%M-0600')

    # Set index
    df.set_index('datetime', inplace=True)
    df.sort_index(inplace=True)

    # Drop unnecessary columns
    # df.drop(['POC', 'Date Local', 'Time Local',
    #          'Date GMT', 'Time GMT', 'Uncertainty', 'Qualifier',
    #          'Method Type', 'Method Code', 'Method Name',
    #          'Date of Last Change'], axis=1, inplace=True)

    # Create extra location field for avoiding duplicate entries when pivoting
    # df['County-Site Code'] = df['County Code'] * 10000 + df['Site Num']

    return df


def print_sites():
    df = read_dataframe()
    df_sites = pd.unique(df[['site', 'lat', 'lon', 'GISDatum', 'elev']].values)
    print(df_sites)


def test_site():
    df = read_dataframe()
    df_current = df[(df['site'] == '840010970003') & (df['parameter'] == 88502.0)]

    ts = pd.Series(df_current['value'])
    print(ts.head())

    plt.figure()
    ts.plot()
    plt.show()


if __name__ == "__main__":
    df = read_dataframe()
    df_parameters = pd.unique(df[['parameter']].values)

    for param in df_parameters:
        print('Processing Parameter: {}'.format(param[0]))

        df_current = df[(df['parameter'] == param[0])]
        if df_current.size == 0:
            print("No rows!")
            continue

        try:
            pivot = df_current.pivot(columns='site',
                                     values='value')

            fig = plt.gcf()
            fig.set_size_inches(12, 7)
            ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=2)
            pivot.plot(ax=ax1)
            ax1.set_title(param)

            ax2 = plt.subplot2grid((1, 3), (0, 2))
            pivot.plot.hist(alpha=0.5, ax=ax2)

            fig.tight_layout()
            plt.show()
        except Exception as e:
            print(e)
