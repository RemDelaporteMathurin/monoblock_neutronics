from run_neutronics import my_source

import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sigma):
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

def muir_spectrum(source):
    E = source.energy
    mu = E.e0
    sigma = np.sqrt(4.*E.e0*E.kt/E.m_rat)
    return mu, sigma

def plot_muir_spectrum(source, n_samples, **kwargs):
    mu, sigma = muir_spectrum(source)
    energy = np.random.normal(mu, sigma, n_samples)
    energy *= 1e-6  # eV to MeV

    return plt.hist(energy, **kwargs)

n_samples = 100*50
n, bins, patches = plot_muir_spectrum(my_source, n_samples, bins="auto", alpha=0.8, density=True, edgecolor='white')
mu, sigma = muir_spectrum(my_source)
x = np.linspace(11e6, 17e6, num=1000)

plt.plot(x*1e-6, gaussian(x*1e-6, mu*1e-6, sigma*1e-6))
x_annotation = (mu+sigma)*1e-6
y_annotation = gaussian(x_annotation, mu*1e-6, sigma*1e-6)
x_annotation += 0.5

plt.annotate("$\mu = E_0 = $ {:.0f} MeV".format(mu*1e-6), (x_annotation, y_annotation))
plt.annotate(r"$\sigma = \sqrt{\frac{4 E_0 T_i}{M_\mathrm{reactants}}}$", (x_annotation, y_annotation*0.8))


plt.xlabel("Neutron energy (MeV)")
plt.ylabel("Probability")

plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)

plt.savefig("out.png")
plt.savefig("out.svg")
plt.savefig("out.pdf")
