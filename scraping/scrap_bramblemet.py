import datetime
import urllib.request
import re
from bs4 import BeautifulSoup
from pymongo import MongoClient

#####################################################33
### Data scraping
# url="http://www.bramblemet.co.uk/archive/2017/June/CSV/Bra13Jun2017.csv"
# urllib.request.urlretrieve(url, 'data/Bramblemet/scraping/Bra13Jun2017.csv')


# Urls
URLS = ["http://www.bramblemet.co.uk/archive/{}/{}/CSV/Bra{}.csv",
        "http://www.cambermet.co.uk/archive/{}/{}/CSV/Cam{}.csv",
        "http://www.chimet.co.uk/archive/{}/{}/CSV/Chi{}.csv",
        "http://www.sotonmet.co.uk/archive/{}/{}/CSV/Sot{}.csv"]

# Output folder
OUTPUT_FOLDERS = ["Bramblemet", "Cambermet", "Chimet", "Sotonmet"]

# Scraping
option = 3;
date = datetime.date(2017, 1, 1)
while date < datetime.date(2017, 6, 14):
    print("Scraping date: {}".format(date))

    try:
        url = URLS[option] \
            .format(date.strftime("%Y"),
                    date.strftime("%B"),
                    date.strftime("%d%b%Y"))
        print(url)

        filename = "{}.csv".format(date.strftime("%Y-%m-%d"))
        urllib.request.urlretrieve(url, 'data/Bramblemet/scraping/{}/{}'.format(
            OUTPUT_FOLDERS[option], filename))
    except urllib.request.HTTPError as e:
        print(e)

    date += datetime.timedelta(days=1)

# Join files (Using bash script)
