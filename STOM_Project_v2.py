# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:11:15 2020

@author: dd719
"""
#%%
#the firse cell is dedicated to importing all relevant packages, if a package is needed later down the line, add it to the list here
import STOM_higgs_tools#this line imports all the funcyions contained in the file avaidible on bb, for this to work, make sure that file is visible in the files window in your ide
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
#new comment: every time you import STOM_higgs_tools you get a different data set, make sure you make your code general and don't use any values from your specific data set as they will change
#sometimes the data you get looks really bad and you seemingly get a second peak, if that happens to you, re-run every cell in order until you get a graph that looks reasonable
#%%
print('Part 1')
data_array=STOM_higgs_tools.generate_data()
#print(data_array)#there will likeley be a lot of prints for my peace of mind to chech what im working with at each point
#the sheet recommends a histogram function from matplotlib, but i will write my own as it will be easier to extract data later down the line, and can be made to look like their figure
n_bins=30#self-explanatory, this determines the number of bins to be made
bin_boundary_array=np.linspace(104,155,(n_bins+1))#makes an array of consecutive bin boundaries
#print(bin_boundary_array)
bin_width=(155-104)/n_bins
#print(bin_width)
bins = [[] for _ in range(n_bins)]#makes the required number of empty arrays to act as bins
#print(bins)
counter1=0
for i in bin_boundary_array:#this is a lot of comparisons made and may take a couple seconds to process
    if counter1 == n_bins:
        break# this is vital as an error occurs when the loop looks for a nonexistant index of 31 in bin_boundary_array
    for i in data_array:
        if i >= bin_boundary_array[counter1]:
            if i < bin_boundary_array[counter1+1]:#notice the >= and < ensure that edge values get added exactly once
                bins[counter1].append(i)
    counter1=counter1+1
#print(len(bins))#yup, its equal to n_bins
#################################
#this cell will create our x coordinates, the bin centres, and our y coordinates, the number of values in each bin
counter2=0
bin_centre_array=[]
for i in bin_boundary_array:#we aren't actually operating on i, we just want the correct number of operations
    if counter2 == n_bins:
        break#as above
    bin_centre_array.append((bin_boundary_array[counter2]+bin_boundary_array[counter2+1])/2)
    counter2=counter2+1
#print(len(bin_centre_array))#nice, still equal to n_bins
entry_number_array=[]
for i in bins:
    entry_number_array.append(len(i))
#print(len(entry_number_array))#nice, still equal to n_bins
#################################
#so they want errors it turns out, what these are is summerised in the project diary word document, this cell will find them
#first, the slightly easier as it turns out, height error:
height_error=[]
for i in bins:
    height_error.append(np.sqrt(len(i)-((len(i))**2)/len(data_array)))
#now, the error in the energy
energy_error=[]
for i in bins:
    current_bin=np.array(i)
    sum_of_squares=sum(current_bin*current_bin)
    varience=(sum_of_squares/len(i))-(sum(current_bin)/len(current_bin))**2
    energy_error.append(np.sqrt(varience))
#print(energy_error)
#################################
#now for plotting, and some visual stuff at the end:
#entry_number_array=np.array(entry_number_array)/bin_width
#height_error=np.array(height_error)/bin_width     #these would give freq. density

plt.errorbar(bin_centre_array,entry_number_array,xerr=energy_error,yerr=height_error,fmt='o', mew=1, ms=1, capsize=3)
plt.hist(data_array, range = [104, 155], bins = 30, color = "#4dc6ff")
plt.plot(bin_centre_array,entry_number_array,'.',c='#9300FF')

#the above line plots our points, everything else is vanity
plt.grid()
plt.xlabel("$m_{\gamma \gamma } (GeV)$") 
plt.ylabel("Number of Entries") 
plt.title("Photon Energy Histogram") 
plt.legend(["Data",'Error','Data'], loc=1 ) 
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
plt.grid(b=True, which='major', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.savefig('part 1 grpph')
plt.show()

#################################################################################
print('Part 2')

from scipy.optimize import curve_fit

heights, edges, _ = plt.hist(data_array, range = [104, 155], bins = 30)
plt.clf()

data_array.sort()
for i in range(len(data_array)):
    if data_array[i] > 120:
        index_120 = i - 1
        break

bckgrnd = data_array[: index_120]

bin_width = edges[1] - edges[0]
range_ = 120

n_bins = int(range_ // bin_width)

y, x, _ = plt.hist(bckgrnd, bins = n_bins, color='#4dc6ff')
x = x[:-1]
x += bin_width / 2

def to_fit(x, A, lm):
    return A * np.exp(-x / lm)

p0 = [160000, 30]

fit_params, _ = curve_fit(to_fit, x, y, p0 = p0)

fit = lambda x: to_fit(x, *fit_params)
plt.grid(b=True, which='major', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.plot(x, list(map(fit, x)), color = "#e100ff")
plt.title("Background Fit")
plt.xlabel("$m_{\gamma \gamma } (GeV)$") 
plt.ylabel("Frequancy Density") 
plt.legend(['Fit',"Data"], loc=1 ) 
plt.rcParams["figure.figsize"] = fig_size
plt.savefig('part 2 grpph')
plt.show()

#################################################################################
print("Part 3")

from scipy.stats.distributions import chi2

chi_squared_bk = STOM_higgs_tools.get_B_chi(bckgrnd, (0, 120), n_bins, fit_params[0], fit_params[1])
print(f"\tThe reduced chi squared for the background only is {chi_squared_bk}")
print(f"\tThis has a corresponding p value of {chi2.sf(chi_squared_bk, n_bins - 2)}")

#################################################################################

print("Part 4")

for i in range(len(data_array)):
    if data_array[i] > 104:
        index_104 = i
        break

for i in range(len(data_array)):
    if data_array[i] > 155:
        index_155 = i
        break

signal = data_array[index_104: index_155]

chi_squared_signal = STOM_higgs_tools.get_B_chi(signal, (104, 155), 30, fit_params[0], fit_params[1])
print(f"\tThe reduced chi squared for the singal is {chi_squared_signal}")
print(f"\tThis has a corresponding p value of {chi2.sf(chi_squared_signal, 30 - 2)}")
print("\tThis is significant enough to reject the background only (null) hypothesis")
print("\tBy trial and error, the signal amplitude that yields a p-value of 0.05 is about 220")

plt.plot(edges[:-1] + bin_width // 2, list(map(fit, edges[:-1] + bin_width // 2)), color = "#e100ff")
plt.hist(data_array, range = [104, 155], bins = 30, color = "#4dc6ff")

plt.errorbar(bin_centre_array,entry_number_array,xerr=energy_error,yerr=height_error,fmt='o', mew=1, ms=1, capsize=3,c='#0022ff')
plt.plot(bin_centre_array,entry_number_array,'.',c='#ff002f')
#the above line plots our points, everything else is vanity
plt.grid()
plt.xlabel("$m_{\gamma \gamma } (GeV)$") 
plt.ylabel("Number of Entries") 
plt.title("Photon Energy Histogram") 
plt.legend(["Fit",'Data','Data','Error'], loc=1 ) 
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
plt.grid(b=True, which='major', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.savefig('photon energy histogram plus fit')
plt.savefig('part 4 grpph')
plt.show()

############################################################################################
print("Part 5")
def get_SB_fit(x, A, lamb, mu, sig, signal_amp):
        return A*np.exp(-x/lamb) + STOM_higgs_tools.signal_gaus(x, mu, sig, signal_amp)
    
p0_sb=[160000,30,125,1.5,700]
sb_fit_params,_= curve_fit(get_SB_fit,edges[:-1] + bin_width // 2,heights,p0=p0_sb)
fit_sb=lambda x: get_SB_fit(x, *sb_fit_params)
plt.plot(edges[:-1] + bin_width // 2, list(map(fit_sb, edges[:-1] + bin_width // 2)),color='#e100ff')
plt.hist(data_array, range = [104, 155], bins = 30, color = "#4dc6ff")
plt.title("Background+Signal fit")



plt.errorbar(bin_centre_array,entry_number_array,xerr=energy_error,yerr=height_error,fmt='o', mew=1, ms=1, capsize=3,c='#0022ff')
plt.plot(bin_centre_array,entry_number_array,'.',c='#ff002f')
#the above line plots our points, everything else is vanity
plt.grid()
plt.xlabel("$m_{\gamma \gamma } (GeV)$") 
plt.ylabel("Number of Entries") 
plt.title("Photon Energy Histogram") 
plt.legend(["Fit",'Data','Data','Error'], loc=1 ) 
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
plt.grid(b=True, which='major', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.savefig('photon energy histogram plus fit')
plt.savefig('part 5 grpph')
plt.show()


def get_SB_chi(vals, mass_range, nbins, A, lamb, mu, sig, signal_amp):
    bin_heights, bin_edges = np.histogram(vals, range = mass_range, bins = nbins)
    half_bin_width = 0.5*(bin_edges[1] - bin_edges[0])
    ys_expected = get_SB_fit(bin_edges + half_bin_width, A, lamb, mu, sig, signal_amp)
    chi = 0

    # Loop over bins - all of them for now. 
    for i in range( len(bin_heights) ):
        chi_nominator = (bin_heights[i] - ys_expected[i])**2
        chi_denominator = ys_expected[i]
        chi += chi_nominator / chi_denominator
    
    return chi #/float(nbins-2) # B has 2 parameters.

chi_squared_signal_background = get_SB_chi(signal, (104, 155), 30, sb_fit_params[0], sb_fit_params[1],sb_fit_params[2],sb_fit_params[3],sb_fit_params[4])
print(f"\tThe reduced chi squared for the background only is {chi_squared_signal_background}")
print(f"\tThis has a corresponding p value of {chi2.sf(chi_squared_signal_background, n_bins - 2)}")



