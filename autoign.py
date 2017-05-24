import sys
import numpy as np
from cantera import *
import csv

def max_dt(vec, nt):
        nmax = 0
        dmax = 0
        for n in range(nt-2):
                dvec = (vec[n+1] - vec[n])
                if (dvec > dmax):
                        nmax = n
                        dmax = dvec
        return nmax
        
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def find_auto_ign(gas): 
        #steps
        nt = 3000
        dt = 1.e-3 / nt#s
        #Storage space
        #mfrac = {} #tu beda sklady masowe gazow
        #for i in range(0, nt):
        #        mfrac[i] = np.zeros(gas.n_species,'d') #wektor wektorow
        tim = np.zeros(nt,'d') #wektor o wymiarze nt
        temp = np.zeros(nt,'d')
        dtemp = np.zeros(nt-1,'d')
        pressure = np.zeros(nt, 'd')
        #################################################################
        #Create the batch reactor
        r = IdealGasReactor(gas)
        # Now create a reactor network consisting of the single batch reactor
        sim = ReactorNet([r])
        #Run the simulation
        # Initial simulation time
        time = 0.0
        #Loop for nt time steps of dt seconds.

        for n in range(nt):
                time += dt
                sim.advance(time)
                tim[n] = time
                temp[n] = r.T
                #mfrac[n][:] = r.thermo.Y
                pressure[n] = r.thermo.P #Pa
                #r.thermo.Y - sklad masowy, .X - sklad molowy
                #http://www.cantera.org/docs/sphinx/html/cython/thermo.html
                
        #################################################################
        # Catch the autoignition timing
        #################################################################
        
        Dtmax = max_dt(temp, nt)
        Autoignition = tim[Dtmax + 1]

        pressure_smooth = smooth(pressure, 15)
        Dpmax = max_dt(pressure, nt)
        dpdt_max = (pressure_smooth[Dpmax + 1] - pressure_smooth[Dpmax])/(dt) 
        K = dpdt_max / 1e6 * 1#m3; K = MPa*m/s

        Autoignition = Autoignition*1000 #ms
        FinalTemp = temp[-1]

        p_max = max(pressure)
        
        
        #print "species_names: " + ', '.join(r.thermo.species_names)
        #index = r.thermo.species_index('C2H2') #indeks w wektorze mfrac
        #to_plot = np.zeros(nt,'d')
        #for i in range(0, nt):
        #        to_plot[i] = mfrac[i][index]
        #plot(tim, to_plot)
        #show()
        
        #print(r.thermo.report())
 
        return [Autoignition, FinalTemp, K, p_max]

#main

print "Computing:"

data = []
rg = 20
for i in range(rg + 1):
        gas = Solution('gri30.cti')
        
        theta = 1.*i / rg
        
        #N2_mole == 0: 7.52 * O2pzero - O2_mole + C2H2_mole * 2.5 == 0
        # 7.52 *  C2H2_mole * 2.5 / 2 - C2H2_mole * 2.5 / min_phi  + C2H2_mole * 2.5 == 0
        # 7.52 *  2.5 / 2 - 2.5 / min_phi  + 2.5 == 0
        # 2.5/ min_phi == 7.52 * 2.5/2 + 2.5
        # min_phi = 2.5/(7.52 * 2.5/2 + 2.5)
        min_phi = 1 / (7.52 / 2 + 1)
        
        phi = 1 - theta * ( 1 - min_phi)
        C2H2_mole = 2
        
        #yCmHn + y(m + n/4)O2 -> ...
        #m = 2; n = 2; y = C2H2_mole; 
        #O2_mole = C2H2_mole * (m +n /4) / phi

        O2_mole = C2H2_mole * 2.5 / phi
        
        O2pzero = C2H2_mole * 2.5 / 2
        
        N2_mole = 7.52 * O2pzero - O2_mole + C2H2_mole * 2.5 
        H2O_mole = 0
        
        
        
        s = 'C2H2:' + str(C2H2_mole) + ',O2:' + str(O2_mole) + ',N2:' + str(N2_mole) + ',H2O:' + str(H2O_mole)
        print(s)
        
        gas.TPX = 1200, one_atm, s
        
        data.append([phi, C2H2_mole, O2_mole, N2_mole, H2O_mole] + find_auto_ign(gas))

print "done"


csv_file = 'data.csv'
with open(csv_file, 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['phi', 'C2H2_mole', 'O2_mole', 'N2_mole', 'H2O_mole', 'Auto ignition time [ms]','Final Temperature [K]','K [MPa*m/s]', 'Max pressure Pa'])
        for d in data:
                writer.writerow(d)
        print 'output written to '+csv_file
