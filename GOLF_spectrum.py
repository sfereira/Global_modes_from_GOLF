
#-----------------------------------------------------------------------------------------------------------#
#-- Libraries required --#                                                                                  #
#                                                                                                           #
from matplotlib.backend_bases import MouseButton                                                            #
import matplotlib.pyplot as plt                                                                             #
import numpy as np                                                                                          #
import matplotlib                                                                                           #
from matplotlib.lines import Line2D                                                                         #
from astropy.io import fits                                                                                 #
import math                                                                                                 #
from scipy import fftpack                                                                                   #
from scipy import optimize                                                                                  #
from scipy import signal                                                                                    #
from scipy.ndimage import gaussian_filter                                                                   #
from scipy.signal import savgol_filter                                                                      #
import matplotlib.patches as mpatches                                                                       #
#                                                                                                           #
#-----------------------------------------------------------------------------------------------------------#
#-- Reading the input file --#                                                                              #
                                                                                                            #
hdul = fits.open('GOLF_22y_MEAN.fits')                                                                      #
hdul.info()                                                                                                 #
v = hdul[0].data                                                                                            #
#-----------------------------------------------------------------------------------------------------------#
#-- The number of data points --#                                                                           #
                                                                                                            #
n = len(v)                                                                                                  #
print('The number of data points (n) = ', n )                                                               #
#-----------------------------------------------------------------------------------------------------------#
#-- Fast fourier transform of velocity of solar oscillation --#                                             #
                                                                                                            #
a= np.fft.fft(v)                                                                                            #
#-----------------------------------------------------------------------------------------------------------#
#-- Plot a sample of the velocity data --#                                                                  #
                                                                                                            #
print('\t ***--- Plot a sample of the velocity data ---*** \t ')                                            #
fig, ax = plt.subplots()                                                                                    #
ax.set_title('Plot of 1st 100 velocity data')                                                               #
ax.plot(v[0:100])                                                                                           #
plt.show()                                                                                                  #
                                                                                                            #
ax.set_title('Plot of 1st 100 amplitudes')                                                                  #
ax.plot(a[0:100])                                                                                           #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Power Calculation                                                                                         #
                                                                                                            #
power = ( abs( a )**2 )                                                                                     #
power = power * 1e4                                                                                         #
                                                                                                            #
#-- Plot a sample of the power data --#                                                                     #
plt.title('Plot of 1st 100 power data')                                                                     #
plt.plot(power[0:100])                                                                                      #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
#-- Number of independent frequency points --#                                                              #
                                                                                                            #
n2 = round(n / 2)                                                                                           #
print('Number of independent frequency points =', n2)                                                       #
#-----------------------------------------------------------------------------------------------------------#
#-- Maximum frequency numax = 1. / (2*deltat)*100 (mHz) --#                                                 #
                                                                                                            #
deltat = 20.                                                                                                #
numax = 1. / (2. * deltat) * 1000                                                                           #
print('Maximum frequency numax [mHz] = ',numax)                                                             #
#-----------------------------------------------------------------------------------------------------------#
#- Define frequency points -#                                                                               #
                                                                                                            #
nu = np.arange(n2)/(n2 - 1.)*numax                                                                          #
#-----------------------------------------------------------------------------------------------------------#
#-- l=0, n=14 line --#                                                                                      #
#-- Plot of Power vs frequency --#                                                                          #
                                                                                                            #
fig, ax = plt.subplots()                                                                                    #
ax.set_xlim(2.08,2.10)                                                                                      #
ax.set_ylim(0.000,3e13)                                                                                     #
ax.plot(nu,power[0:n2])                                                                                     #
                                                                                                            #
coords=[]                                                                                                   #
q = 0                                                                                                       #
                                                                                                            #
#-- Function to print mouse click event coordinates --#                                                     #
                                                                                                            #
def onclick(event):                                                                                         #
#    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %                                          #
#          ('double' if event.dblclick else 'single', event.button,                                         #
#           event.x, event.y, event.xdata, event.ydata))                                                    #
    global x1, y1, q                                                                                        #
    x1 = event.xdata                                                                                        #
    y1 = event.ydata                                                                                        #
    if q == 0 :                                                                                             #
        print('select the line frequency range, click to choose left bound and then the right one')         #
        q +=2                                                                                               #
    if q == 2 :                                                                                             #
        print('select the line width, click to choose left wing and then the right one')                    #
        q +=2                                                                                               #
    if q == 4 :                                                                                             #
        print('select the initial mode frequency and amplitude')                                            #
        q+=1                                                                                                #
                                                                                                            #
    print('x = %f, y = %f' %(event.xdata, event.ydata))                                                     #
                                                                                                            #
    global coords                                                                                           #
    coords.append((event.xdata, event.ydata))                                                               #
                                                                                                            #
    if len(coords) == 5:                                                                                    #
        fig.canvas.mpl_disconnect(cid)                                                                      #
                                                                                                            #
    return coords                                                                                           #
                                                                                                            #
cid = fig.canvas.mpl_connect('button_press_event', onclick)                                                 #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
#-- Data after selecting new range --#                                                                      #
                                                                                                            #
nul = coords[0][0]                                                                                          #
nur = coords[1][0]                                                                                          #
                                                                                                            #
wl = coords[2][0]                                                                                           #
wr = coords[3][0]                                                                                           #
                                                                                                            #
nuj0 = coords[4][0]                                                                                         #
aj0 = coords[4][1]                                                                                          #
                                                                                                            #
wj0 = wr - wl                                                                                               #
                                                                                                            #
c0=0                                                                                                        #
c1=0                                                                                                        #
                                                                                                            #
x0 = [aj0,nuj0,wj0,c0,c1]                                                                                   #
print('\n')                                                                                                 #
print('\t Initial values \t')                                                                               #
print('\n')                                                                                                 #
print('[ amplitude, mode frequency, linewidth, noise_coefficients ]', x0)                                   #
print('\n')                                                                                                 #
                                                                                                            #
il = np.where(min(abs(nu-nul)) == abs(nu-nul))                                                              #
ir = np.where(min(abs(nu-nur)) == abs(nu-nur))                                                              #
                                                                                                            #
print(il,ir)                                                                                                #
                                                                                                            #
nu1 = nu[il[0][0]:ir[0][0]]                                                                                 #
                                                                                                            #
pow1 = power[il[0][0]:ir[0][0]]                                                                             #
                                                                                                            #
n1 = len(nu1)                                                                                               #
print('n1=',n1)                                                                                             #
                                                                                                            #
x0= np.array([aj0,nuj0,wj0,c0,c1])                                                                          #
                                                                                                            #
ftol=1e-8                                                                                                   #
                                                                                                            #
xi = np.zeros([5,5])                                                                                        #
                                                                                                            #
for i in range(0,4):                                                                                        #
    xi[i,i] = 1.                                                                                            #
                                                                                                            #
x00 = x0                                                                                                    #
                                                                                                            #
#-----------------------------------------------------------------------------------------------------------#
# function used in the minimizing in the Powell method.                                                     #
# Lorentzian line profile plus linear noise.                                                                #
# aj - amplitude                                                                                            #
# wj - line width                                                                                           #
# nuj - mode frequency                                                                                      #
# c0, c1 - noise coefficients                                                                               #
# s - log of maximum likelihood function                                                                    #
#-----------------------------------------------------------------------------------------------------------#
def sfunc(x):                                                                                               #
    aj = x[0]                                                                                               #
    nuj=x[1]                                                                                                #
    wj=x[2]                                                                                                 #
    c0=x[3]                                                                                                 #
    c1=x[4]                                                                                                 #
    wj2=(wj/2.)**2                                                                                          #
    s=0                                                                                                     #
    for i in range(n1-1):                                                                                   #
        mi = aj * wj2  / ( (nu1[i] - nuj)**2 + wj2 ) + c0 + c1 * nu1[i]                                     #
        if mi > 0:                                                                                          #
            s = s + math.log(mi) + pow1[i]/mi                                                               #
                                                                                                            #
    return s                                                                                                #
#-----------------------------------------------------------------------------------------------------------#
# Minimizing function                                                                                       #
                                                                                                            #
minimum = optimize.fmin_powell(sfunc, x0, ftol=ftol, maxiter=500)                                           #
                                                                                                            #
x=minimum                                                                                                   #
aj=x[0]                                                                                                     #
nuj=x[1]                                                                                                    #
wj=x[2]                                                                                                     #
c0j=x[3]                                                                                                    #
c1j=x[4]                                                                                                    #
wj2=(wj/2.)**2                                                                                              #
                                                                                                            #
plt.title('GOLF Power spectrum 1996-2018')                                                                  #
plt.xlabel(r'$\nu$  (mHz)')                                                                                 #
plt.ylabel(r'Power')                                                                                        #
plt.plot(nu1,pow1)                                                                                          #
plt.plot(nu1, aj*(wj2/((nu1-nuj)**2 + wj2 )) + c0j + c1j * nu1,color='red')                                 #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Smoothening (By convolve) the results to reduce noise and see better results.                             #
                                                                                                            #
win = signal.windows.hann(200,sym=False)                                                                    #
w = win/2                                                                                                   #
                                                                                                            #
filtered = signal.convolve(pow1, w, mode='same',method='fft')/((sum(win)/10))                               #
                                                                                                            #
plt.title('GOLF Power spectrum (smoothed) - Convolution')                                                   #
plt.xlabel(r'$\nu$  (mHz)')                                                                                 #
plt.ylabel(r'Power')                                                                                        #
plt.plot(nu1,filtered*1.1)                                                                                  #
plt.plot(nu1, aj*wj2/((nu1-nuj)**2 + wj2 ) + c0j + c1j * nu1)                                               #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Smoothening (By gaussian) the results to reduce noise and see better results.                             #
                                                                                                            #
gf = gaussian_filter(pow1, sigma=5)*sum(win)/19                                                             #
                                                                                                            #
fig, ax = plt.subplots()                                                                                    #
ax.set_title('GOLF Power spectrum (smoothed) - gaussian_filter ')                                           #
ax.set_xlabel(r'$\nu$  (mHz)')                                                                              #
ax.set_ylabel(r'Power')                                                                                     #
ax.plot(nu1,gf*1.1)                                                                                         #
ax.plot(nu1, aj*wj2/((nu1-nuj)**2 + wj2 ) + c0j + c1j * nu1)                                                #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Smoothening (By Savgol filter) the results to reduce noise and see better results.                        # 
                                                                                                            #
sf = savgol_filter(pow1,221,1,mode='wrap')*sum(win)/20                                                      #
                                                                                                            #
fig, ax = plt.subplots()                                                                                    #
ax.set_title('GOLF Power spectrum (smoothed) - savgol_filter ')                                             #
ax.set_xlabel(r'$\nu$  (mHz)')                                                                              #
ax.set_ylabel(r'Power')                                                                                     #
ax.plot(nu1,sf)                                                                                             #
ax.plot(nu1, aj*wj2/((nu1-nuj)**2 + wj2 ) + c0j + c1j * nu1)                                                #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Plotting and comparing the results from all the smoothening.                                              #
                                                                                                            #
fig, ax = plt.subplots()                                                                                    #
ax.plot(nu1,filtered,color='red')                                                                           #
ax.plot(nu1,gf,color='blue')                                                                                #
ax.plot(nu1,sf,color='green')                                                                               #
ax.plot(nu1,pow1,color='black')                                                                             #
ax.plot(nu1, aj*wj2/((nu1-nuj)**2 + wj2 ) + c0j + c1j * nu1,color='purple')                                 #
                                                                                                            #
red_patch = mpatches.Patch(color='red', label='Convolution')                                                #
blue_patch = mpatches.Patch(color='blue', label='gaussian_filter')                                          #
green_patch = mpatches.Patch(color='green', label='savgol_filter')                                          #
                                                                                                            #
plt.legend(loc='upper right',fontsize='large',handles=[red_patch,blue_patch,green_patch])                   #
plt.show()                                                                                                  #
#-----------------------------------------------------------------------------------------------------------#
