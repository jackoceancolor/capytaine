# script to plot RM3 data
import matplotlib.pyplot as plt

plt.figure()
plt.plot(rm3_results.omega, rm3_results.wavelength,
     rm3_results.omega, (9.81/(rm3_results.omega**2))**1)
plt.ylabel('wavenumber')
plt.xlabel('omega')
plt.legend({'wavenumber', 'mycalc'})
plt.ylim(0, 100)
plt.show()

# input("hit[enter] to end.")
# plt.close()
