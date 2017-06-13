bramblemet <- read.csv('data/Bramblemet/BraTest.csv',header=T)
#bramblemet <- read.csv('data/Bramblemet/Sot06Jun2017.csv',header=T)
bramblemet$Date <- paste(bramblemet$Date,bramblemet$Time)
bramblemet$Date <- strptime(bramblemet$Date,format="%d/%m/%Y %H:%M")
bramblemet <- bramblemet[,-2]
plot(bramblemet)

#bramblemet = bramblemet[bramblemet$Date>='2017-06-06' & bramblemet$Date<'2017-06-07',]
plot(bramblemet$WSPD, type='l')

plot(bramblemet$WSPD, type='l')

plot(bramblemet$WD, type='l')
plot(bramblemet$GST, type='l')
plot(bramblemet$ATMP, type='l')
plot(bramblemet$WTMP, type='l')
plot(bramblemet$BARO, type='l')
plot(bramblemet$DEPTH, type='l')
plot(bramblemet$VIS, type='l')
