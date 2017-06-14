import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_dataframe():
    # Read input file
    # Rain hourly
    df = pd.read_csv('data/midas/midas_rainhrly_201701-201712.txt', header=None,
                     parse_dates=[0])

    # Pre-processing
    #df[0] = pd.to_datetime(df[0], format='%Y-%m-%d %H:%M')
    df[2] = [x.strip() for x in df[2]]
    df[5] = [x.strip() for x in df[5]]

    def my_parse_float(value):
        try:
            return float(value)
        except ValueError:
            return np.NAN

    df[8] = df[8].apply(lambda x: my_parse_float(x))

    # Drop unnecessary columns
    #df.drop([1, 3, ], axis=1, inplace=True)

    # Create extra location field for avoiding duplicate entries when pivoting
    #df['County-Site Code'] = df['County Code'] * 10000 + df['Site Num']

    return df


if __name__ == "__main__":
    df = read_dataframe()
    #print(df[0].head())
    #exit(0)

    df_current = df[(df[1] == 1586) & (df[2] == 'RAIN') & (df[3] == 1) & (df[4] == 1) & (df[5] == 'SREW')]
    df_current.set_index(0, inplace=True)
    if df_current.size == 0:
        print("No rows!")

    #print(df_current.head())
    #exit(0)

    ts = pd.Series(df_current[8])
    print(ts.head())
    #exit()

    ts.plot()
    #plt.title(state[1])
    #plt.ylabel("{} [{}]".format(param[1],param[2]))
    plt.show()
