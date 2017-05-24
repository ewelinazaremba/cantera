import csv
from matplotlib.pylab import *

csv_file = 'data.csv'

#data
data = [ [] ]
frac = []
autoign = []
temp = []
K = []
press = []

#read and parse
with open(csv_file, 'r') as infile:
        reader = csv.reader(infile)
        i = 0
        for row in reader:
                i = i + 1
                if i == 1:
                        for n in range(len(data), len(row)):
                                data.append([])
                        continue
                
                for j, elem in enumerate(row):
                        val = float(elem)
                        data[j].append(val)
                        
frac = data[0]
autoign = data[5]
temp = data[6]
K = data[7]
press = data[8]
for i in range(0, len(press)):
        press[i] /= 1000
                

f, axarr = plt.subplots(2, 2)
axarr[0][0].plot(frac, autoign, '.')
axarr[0][0].set_title("Autoignition [ms]")
axarr[0][0].set_xlabel("phi")
axarr[0][1].plot(frac, temp, '.')
axarr[0][1].set_title("Final temp [K]")
axarr[0][1].set_xlabel("phi")
axarr[1][0].plot(frac, K, '.')
axarr[1][0].set_title("K [MPa*m/s]")
axarr[1][0].set_xlabel("phi")
axarr[1][1].plot(frac, press, '.')
axarr[1][1].set_title("Max pressure [kPa]")
axarr[1][1].set_xlabel("phi")

axarr[0][0].grid()
axarr[0][1].grid()
axarr[1][0].grid()
axarr[1][1].grid()

f.subplots_adjust(hspace=0.5)

#plot(Ti2,Autoignition_cas, '^', color = 'orange')
#xlabel(r'Temp [1000/K]', fontsize=20)
#ylabel("Autoignition [ms]")
#title(r'Autoignition of $CH_{4}$ + Air mixture at $\Phi$ = 1, and P = 1 bar',
#fontsize=22,horizontalalignment='center')
        #
        

#        to_plot = np.zeros(nt,'d')
#        for i in range(0, nt):
#                to_plot[i] = mfrac[i][index]
        
        
#        plot(tim, to_plot)
        #plot(tim, pressure)       
        
        #axis([0, nt*dt, 0, 1])

show()

                
