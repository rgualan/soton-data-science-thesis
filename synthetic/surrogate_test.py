import numpy
import pandas
import matplotlib.pyplot as plt

numpy.random.seed(0)
ts = numpy.random.normal(0, 1, 1000)

plt.figure()
pandas.tools.plotting.autocorrelation_plot(ts)
plt.ylim([-0.1,0.1])
plt.title('Autocorrelation function of random time series')

ts_fourier  = numpy.fft.rfft(ts)
random_phases = numpy.exp(numpy.random.uniform(0,numpy.pi,int(len(ts)/2)+1)*1.0j)
ts_fourier_new = ts_fourier*random_phases
new_ts = numpy.fft.irfft(ts_fourier_new)

plt.figure()
pandas.tools.plotting.autocorrelation_plot(new_ts)
plt.ylim([-0.1,0.1])
plt.title('Autocorrelation function of surrogate time series')


plt.figure()
plt.plot(ts)
plt.plot(new_ts)

plt.figure()
plt.hist(ts)
plt.hist(new_ts)
plt.show()
