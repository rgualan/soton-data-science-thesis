import pandas as pd
import matplotlib.pyplot as plt


def read_dataframe():
    # Read input file
    # Criteria Gasses
    df = pd.read_csv('data/epa_hourly/hourly_42101_2017_co.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_44201_2017_ozone.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_42401_2017_so2.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_42602_2017_no2.csv')
    # Particulates
    # df = pd.read_csv('data/epa_hourly/hourly_81102_2017_pm10.csv')
    # Meteorological
    # df = pd.read_csv('data/epa_hourly/hourly_WIND_2017.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_TEMP_2017.csv')

    # df = pd.read_csv('data/epa_hourly/hourly_PRESS_2017.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_RH_DP_2017.csv')

    # Parse Date Local
    df['Datetime Local'] = pd.to_datetime(
        df['Date Local'] + ' ' + df['Time Local'],
        format='%Y-%m-%d %H:%M')

    # Drop unnecessary columns
    df.drop(['POC', 'Date Local', 'Time Local',
             'Date GMT', 'Time GMT', 'Uncertainty', 'Qualifier',
             'Method Type', 'Method Code', 'Method Name',
             'Date of Last Change'], axis=1, inplace=True)

    # Create extra location field for avoiding duplicate entries when pivoting
    df['County-Site Code'] = df['County Code'] * 10000 + df['Site Num']

    return df


if __name__ == "__main__":
    df = read_dataframe()
    df_states = pd.unique(df[['State Code', 'State Name']].values)
    df_parameters = pd.unique(df[['Parameter Code', 'Parameter Name', 'Units of Measure']].values)

    for state in df_states:
        print("State: {} # {}".format(state[0], state[1]))

        for param in df_parameters:
            print('Processing Parameter: {} # {} # {}'.format(param[0], param[1], param[2]))

            df_current = df[(df['State Code'] == state[0]) & (df['Parameter Code'] == param[0])]
            if df_current.size == 0:
                print("No rows!")
                continue

            try:
                pivot = df_current.pivot(index='Datetime Local', columns='County-Site Code',
                                         values='Sample Measurement')

                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=2)
                pivot.plot(ax=ax1)
                ax1.set_title(state[1])
                ax1.set_ylabel("{} [{}]".format(param[1],param[2]))

                ax2 = plt.subplot2grid((1, 3), (0, 2))
                pivot.plot.hist(alpha=0.5, ax=ax2)


                plt.show()

            except Exception as e:
                print('Exception: {}'.format(e))
