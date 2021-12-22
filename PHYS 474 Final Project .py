# PHYS 474 Final Project

import numpy as np
from astropy.io import fits
from astropy.io import ascii
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt

####################################
# Data loading and processing
####################################

SDSS_filename = "spec-0266-51602-0001.fits"
hdul = fits.open(SDSS_filename)  # load SDSS data
print(hdul.info())

spec_data = hdul["SPECOBJ"].data
#  u g r i z filters
print(spec_data["SPECTROFLUX"])
data = hdul["COADD"].data

print(data["loglam"])

# Plot of spectral data
plt.figure(figsize=(10, 8))
plt.plot(np.power(10, data["loglam"]), data["flux"], "b", label="Luminous Stars")
plt.xlabel("Wavelength")
plt.ylabel("Flux")
plt.grid()
plt.legend()
plt.title("SDSS Spectral data for object " + SDSS_filename)
plt.tight_layout()
plt.show()

print(spec_data["SPECPRIMARY"])

index_table_data = ascii.read("tmj.dat")  # Load Lick indices data
print(index_table_data)

age = index_table_data["col1"]
Z = index_table_data["col2"]

Hb = index_table_data["col16"]
Mg2 = index_table_data["col19"]
Mgb = index_table_data["col20"]
Fe5270 = index_table_data["col21"]
Fe5335 = index_table_data["col22"]
Fe4531 = index_table_data["col14"]
Fe5015 = index_table_data["col17"]
HdA = index_table_data["col4"]
HgA = index_table_data["col10"]

MgFe = np.sqrt(Mgb * (0.72 * Fe5270 + 0.28 * Fe5335))  # Lick indices described in paper
Mg2Fe = 0.6 * Mg2 + 0.4 * np.log(Fe4531 + Fe5015)  # MgFe and Mg2Fe are sensitive to metallicity
HdA_HgA = HdA + HgA  # Hb and HdA_HgA are sensitive to age

metals = np.array([Z, MgFe, Mg2Fe])
metals = metals.T
metals
metals_valid = metals[~np.isnan(metals).any(axis=1)]  # filter out nand values due to negative MgFe indices
metals_valid.shape

plt.figure(figsize=(10, 8))
plt.scatter(metals_valid[:, 1], metals_valid[:, 0], c="b", label="MgFe")
plt.scatter(metals_valid[:, 2], metals_valid[:, 0], c="r", label="Mg2Fe")
plt.xlabel("Index")
plt.ylabel("Metallicity")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 8))
plt.scatter(Hb, age, c="g", label="Hb")
plt.xlabel("Index")
plt.ylabel("Age")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

print(age)

####################################
# Interpolation and plotting
####################################

# age
x_age = Hb
y_age = HdA_HgA
z_age = age
X_age = np.linspace(min(x_age), max(x_age))
Y_age = np.linspace(min(y_age), max(y_age))
X_age, Y_age = np.meshgrid(X_age, Y_age)  # 2D grid for age interpolation
interp_age = LinearNDInterpolator(list(zip(x_age, y_age)), z_age)  # Interpolate age and age indices from ASCII table
D_age = interp_age(X_age, Y_age)

# metallicity
x_metal = metals_valid[:, 1]
y_metal = metals_valid[:, 2]
z_metal = metals_valid[:, 0]
X_metal = np.linspace(min(x_metal), max(x_metal))
Y_metal = np.linspace(min(y_metal), max(y_metal))
X_metal, Y_metal = np.meshgrid(X_metal, Y_metal)  # 2D grid for metallicity interpolation
interp_metal = LinearNDInterpolator(list(zip(x_metal, y_metal)), z_metal)
D_metal = interp_metal(X_metal, Y_metal)

# age plot
plt.figure(figsize=(10, 8))
plt.pcolormesh(X_age, Y_age, D_age, shading="auto")
# plt.plot(x, y, "ok", label="input point")
plt.legend()
cbar = plt.colorbar()
cbar.set_label("Age (Gyr)")
plt.xlabel("Hb")
plt.ylabel("HdA + HgA")
plt.title("Interpolated age as a function of two Lick indices")
plt.tight_layout()
plt.savefig("age_sensitive_indices.jpeg", dpi=300)
plt.show()

# metallicity plot
plt.figure(figsize=(10, 8))
plt.pcolormesh(X_metal, Y_metal, D_metal, shading='auto')
# plt.plot(x, y, 'ok',  label="input point")
plt.legend()
plt.colorbar(label="Metallicity [Z/H]")
plt.xlabel('MgFe')
plt.ylabel('Mg2Fe')
plt.title("Interpolated metallicity as a function of two Lick indices")
plt.tight_layout()
plt.savefig('metal_sensitive_indices.jpeg', dpi=300)
plt.show()
