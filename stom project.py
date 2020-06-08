print("NOTE: Removed the / float(nbins - 2) from STOM_higgs_tools/get_B_chi, as the reduced chi squared cannot be used for chi squared calculations. \nAccording to arXiv, using the reduced chi squared should be avoided: \n\thttps://arxiv.org/pdf/1012.3754.pdf")

#%% PART 1
############################################################################################
print("Part 1")

from STOM_higgs_tools import *

data = generate_data()

heights, edges, _ = plt.hist(data, range = [104, 155], bins = 30)
plt.show()

#%% PART 2
############################################################################################
print("Part 2")

from scipy.optimize import curve_fit
import numpy as np

data.sort()
for i in range(len(data)):
    if data[i] > 120:
        index_120 = i - 1
        break

bckgrnd = data[: index_120]

bin_width = edges[1] - edges[0]
range_ = 120

n_bins = int(range_ // bin_width)

y, x, _ = plt.hist(bckgrnd, bins = n_bins)
x = x[:-1]
x += bin_width / 2

def to_fit(x, A, lm):
    return A * np.exp(-x / lm)

p0 = [160000, 30]

fit_params, _ = curve_fit(to_fit, x, y, p0 = p0)

fit = lambda x: to_fit(x, *fit_params)

plt.plot(x, list(map(fit, x)))
plt.title("Background fit")
plt.show()

#%% PART 3
############################################################################################
print("Part 3")

from scipy.stats.distributions import chi2

chi_squared_bk = get_B_chi(bckgrnd, (0, 120), n_bins, fit_params[0], fit_params[1])
print(f"\tThe reduced chi squared for the background only is {chi_squared_bk}")
print(f"\tThis has a corresponding p value of {chi2.sf(chi_squared_bk, n_bins - 2)}")

#%% PART 4
############################################################################################
print("Part 4")

for i in range(len(data)):
    if data[i] > 104:
        index_104 = i
        break

for i in range(len(data)):
    if data[i] > 155:
        index_155 = i
        break

signal = data[index_104: index_155]

chi_squared_signal = get_B_chi(signal, (104, 155), 30, fit_params[0], fit_params[1])
print(f"\tThe reduced chi squared for the singal is {chi_squared_signal}")
print(f"\tThis has a corresponding p value of {chi2.sf(chi_squared_signal, 30 - 2)}")
print("\tThis is significant enough to reject the background only (null) hypothesis")
print("\tBy trial and error, the signal amplitude that yields a p-value of 0.05 is about 220")

plt.plot(edges[:-1] + bin_width // 2, list(map(fit, edges[:-1] + bin_width // 2)), color = "red")
plt.hist(data, range = [104, 155], bins = 30, color = "blue")
plt.show()

#%% PART 5
############################################################################################
print("Part 5")
def get_SB_fit(x, A, lamb, mu, sig, signal_amp):
        return A*np.exp(-x/lamb) + signal_gaus(x, mu, sig, signal_amp)
    
p0_sb=[160000,30,125,1.5,700]
sb_fit_params,_= curve_fit(get_SB_fit,edges[:-1] + bin_width // 2,heights,p0=p0_sb)
fit_sb=lambda x: get_SB_fit(x, *sb_fit_params)
plt.plot(edges[:-1] + bin_width // 2, list(map(fit_sb, edges[:-1] + bin_width // 2)),color='red')
plt.hist(data, range = [104, 155], bins = 30, color = "blue")
plt.title("Background+Signal fit")
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



















