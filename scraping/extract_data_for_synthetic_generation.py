import pandas as pd
import matplotlib.pyplot as plt


def read_dataframe():
    # Read input file
    # Criteria Gasses
    # df = pd.read_csv('data/epa_hourly/hourly_42101_2017_co.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_44201_2017_ozone.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_42401_2017_so2.csv')
    # df = pd.read_csv('data/epa_hourly/hourly_42602_2017_no2.csv')
    # Particulates
    df = pd.read_csv('../data/epa_hourly/hourly_81102_2017_pm10.csv')
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

    df_subset = df[(df['State Code'] == 1) & (df['Parameter Code'] == 81102)]
    if df_subset.size == 0:
        print("No rows!")

    pivot = df_subset.pivot(index='Datetime Local', columns='County-Site Code',
                             values='Sample Measurement')

    pivot.drop([730023], axis=1, inplace=True)
    pivot.dropna(inplace=True)
    print(pivot.isnull().sum())

    pivot.to_csv("../data/epa_hourly/alabamatest.csv", index=False, header=False)

    # pivot.plot()
    # plt.show()
