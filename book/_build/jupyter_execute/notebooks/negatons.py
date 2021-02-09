<a href="https://colab.research.google.com/github/AvijeetPrasad/laputas/blob/main/notebooks/negatons.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Gravity of Negative Mass: Some Subtle Aspects and Implications

## Summary
The gravitational interaction has mass of only one sign as the source that is conventionally taken to be positive. The existence of negative mass is not precluded in either Newtonian gravity or general relativity. Negative masses have been invoked in many recent papers to understand the dominance of repulsive dark energy that is accelerating the universe. Here we discuss the various gravitating properties of negative mass particles and their implications for black holes, cosmology, and fundamental physics. Inconsistencies in the earlier works are also discussed. We study the possibility of such negative mass particles accounting for accelerated expansion of the universe. We set constraints on the density of these negative mass particles based on astrophysical observations.

# Import relevant packages and constants
from astropy import units as u
from astropy.constants.si import c, G, m_p, M_sun, hbar,h, k_B
from astropy.cosmology import WMAP9 as cosmo
import math
import matplotlib.pyplot as plt
import numpy as np
#c = speed of light
#m_p =  rest mass of proton
#G =  universal gravatational constant
#M_sun = Solar mass
# Hubble's constant = cosmo.H(0)  
# hbar = Reduced Planck's constant

## 1. Introduction
It is an empirical fact that unlike in electromagnetic interactions where we have both positive and negative electric charges, the gravitational interaction has mass of only one sign as the source. This mass (or energy) is taken to be positive, implying that there are no negative masses in nature. In electromagnetism there is a similar asymmetry in the absence of a monopole which has not been observed, which is not prohibited by theory. Negative mass is quite consistent with gravitational physics [1-3], both Newtonian gravity and general relativity do not preclude its existence. There has not been much discussion about the gravitational effects of negative masses and their other properties. However many recent papers have invoked the possible presence of negative masses to understand cosmic conundrums like the dominance of repulsive dark energy (DE) accelerating the universe [4-6]. 


## 2. Implications of negative mass sources
Negative mass particles with a number density $n$, of mass $m$, moving with velocity $v$ would exert a negative pressure, $P\sim -nmv^2$. It has been pointed out that Einstein himself interpreted the cosmological constant as perhaps due to empty space consisting of gravitating negative masses distributed all over interstellar space [4, 7]. As other fundamental interactions involve charges of both signs as sources of the field, some authors have felt it odd that gravitational charges, conventionally called masses, appear to consist of only positive sign. 
The gravitating properties of negative mass particles (negatons), assuming they exist, would be very interesting and in contrast to the usual positive mass particles. To begin with, the equivalence principle would imply that all test particles whether of positive or negative mass would fall in the field of a positive mass source (like the Earth). This is because we have $-ma=GM(-m)/r^2$ , with -m cancelling out (for negatons both inertial and gravitational masses are negative). However two negatons would gravitationally repel each other so that a collection of such particles would distribute themselves, almost homogenously in the Universe. A negative mass source would repel all particles with both positive and negative mass, whereas positive mass source will attract both. A pair consisting of a positive and negative mass would try to chase each other, but problem is more involved if the particles are moving. In this case it has been shown [2] that the particles move away from each other and the proper separation as measured by an observer commoving with one of the particles increases roughly as $\gamma d$, where $\gamma$ is the Lorentz factor, the separation becoming indefinitely large as $v\rightarrow c$. 
There is symmetry between positive and negative masses in many equations in physics. For instance, relativistic wave equations for particles of any spin are form invariant for particles with both positive and negative mass. This is true for Klein-Gordon, Dirac equations as well as linearized massive spin-2 Pauli-Fierz equation. Dirac equation, well known to have negative energy solutions, interpreted as positrons in negative energy sea. In quantum mechanics it is known that the centre of a wave packet of negative energy solutions describes a uniform motion with velocity in opposite direction to momentum thus behaving like negatons [8]. As far as spin-1/2 particles are concerned, the chirality transformation described by $\psi=\gamma_5 \psi$, converts the Dirac equation for positive mass into one for negative mass. If $\psi (x)$ is the solution for the positive mass case, $\gamma_5 \psi(x)$ is the solution for the negative mass case, the chirality transformation being effectively equivalent to mass inversion, $m\rightarrow -m$ [9]. 
Moreover, negative effective masses occur in many physical contexts like semiconductors [10], Casimir effect [11], Hawking radiation [12] modelled as virtual negatons falling into black hole (BH), excitation of phantom DE as negatons [13], Bose-Einstein condensate [14]. However here we shall mainly be concerned with the gravitating properties of negatons and their implications for BH physics and cosmology. 


## 3. Sign of G for negative mass sources
A common error made in several discussions of the gravity of negative masses, including recent papers [4-6, 15] is to retain the same sign for the gravitational constant G. This would lead to inconsistent results. This is easily seen for instance in considering the gravitational self-energy of mass M with a spatial extent $R$, i.e. 

$E=-\frac{(GM^2)}{R} \quad (1)$ 					

If $M$ is positive, $E$ is the negative gravitational binding energy implying an attractive force. However if $M$ is negative, $M^2$ is positive, but to have a positive binding energy appropriate for a repulsive interaction, $G$ must also change sign. Thus we should have a negative $G$ in discussing the gravitational interaction of negatons. 


### 3.1 Einstein’s field equations with negative mass
The change in sign of $G$ is also to be carried over in discussing the general relativity field equations. Contrary to the inconsistent discussion in recent literature [4, 15] the GR field equations for matter consisting of both positive and negative masses is:

$R_{\mu \nu} - \frac{1}{2} g_{\mu \nu} R = \frac{8\pi G}{c^4}T^+_{\mu \nu} - \frac{8\pi G}{c^4}T^-_{\mu \nu} = \kappa T^+_{\mu\nu} - \kappa T^-_{\mu\nu} \quad (2)$

where $T^+_{\mu \nu}$ and $T^-_{\mu \nu}$ are the energy-momentum tensors for positive and negative masses. The pressure term would be $P=+nmv^2$ or $P=-nmv^2$. The Friedmann equation for the scale factor $a(t)$ would be (for $k=0$):

$\left(\frac{\dot{a}}{a}\right)^2=\frac{8\pi G}{3}\rho_+-\frac{8\pi G}{3}\rho_- \quad (3)$  									

(Equal energy densities of positive and negative masses implies a static Universe)
As is known, inflation scenarios with initial exponential expansion are driven by negative pressure of quantum vacuum. The quantum vacuum has an equation of state for positive mass particles, 

$P=-\rho c^2 \quad (4)$  											
This gives rise to an accelerating universe, as seen from:

$\ddot{a}=-\frac{GM}{a^2} =-\frac{G}{a^2}\left(\frac{4\pi}{3} \left(\rho+\frac{3P}{c^2} \right) a^3 \right) \quad (5)$

For equation of state given by equation (4), $\ddot{a}$ reverses sign. 


###3.2 Inflation with negative mass

Can a quantum vacuum of negatons give rise to an accelerated expansion? Here conservation of energy would imply, 

$P-\rho c^2=0 \Rightarrow P=+\rho c^2 \quad (6)$

However if sign of $G$ is also reversed in equation (5) it is seen that $\ddot{a}$ is positive. This implies that in the initial epoch, the quantum vacuum energy of either positive or negative mass particles or both (if present) can give rise to an inflationary scenario. An effective negative energy (mass) can arise in the early Universe by considering torsion in the early Universe via Einstein-Cartan theory [16, 17, 19] and vortex solutions of this theory could lead to DE [18].

One can discuss the gravitational interactions of both positive and negative mass particles consistently if the sign of $G$ is opposite for negatons. This leads to reasonable results. If we consider a body of mass $M$, with a small particle of mass $\Delta m$ slowly brought from infinity to the surface of $M$, net energy is:

$M+\Delta m - \frac{M(\Delta m) G}{r c^2} = M' \quad (7)$

In all the different cases: $M$, $\Delta m$ positive, $M$, $\Delta m$ negative, or $M$ positive, $\Delta m$ negative, it turns out that Buchdahl's theorem is satisfied, i.e. the maximum energy which can be extracted is either $+Mc^2$ or $-Mc^2$. 

Negatons will have a negative kinetic energy, which in the relativistic case is $-mc^2/\sqrt{1-v^2/c^2}$, but a positive gravitational potential energy. Thus in the Newtonian case:

$-\frac{1}{2} mv^2 + \frac{GMm}{R}=0 \quad (8)$

i.e. negative masses can escape from positive or negative mass sources with escape velocity and have trajectories prescribed by conics. Negatons always have $v\leq c$, like usual positive mass particles but with $m$ replaced by $-m$, a real but negative proper mass. Tachyons in contrast would have imaginary positive proper mass, with $v \geq c$.  


# plot the kinetic energy for different mass
ke = lambda v: -m_p*c**2/np.sqrt(1-(v**2/c**2))
v = np.linspace(.1,.99)*c
plt.plot(v/c,ke(v),'b')
plt.xlabel('v/c')
plt.ylabel('KE')
plt.grid(True)
plt.show()

### 3.3 Negatons as possible candidates for dark matter and dark energy

Gravitationally, dark matter (DM) behaves like baryonic matter, i.e. they have $P=+nm_{DM} v^2$. Therefore DM cannot consist dominantly of negatons. In any case the flux of negatons required to significantly affect masses and evolution of supermassive black holes and even large objects (large capture cross section) like red supergiants is very large. For a billion solar mass black hole, the flux required to significantly alter its mass over the Hubble time ($t_H \approx 4.45\times10^{17} s$) is given as:

$f=\frac{M_{BH}}{(m \times t_H \times A_{BH})}  \approx 2.44\times10^{22} m^{-2} s^{-1} \quad (9)$

where, $A_{BH} \approx 10^{26} m^2$ is the horizon area of the black hole, and the negaton mass, $m$ is taken to be of the order of the proton mass. And similarly we can work out the flux required to alter the evolution of red supergiants as, $f \approx 10^{16} m^{-2} s^{-1}$. These are much higher than the expected DM density of $\approx 2\times 10^{-27} kgm^{-3}$. So negative mass particles cannot be a source of DM. And in any case their energy density being negative, they will not be suitable candidates.  

m_bh = 1.e9 * M_sun 
a_bh = (16 * math.pi * G**2 * m_bh**2) / c**4
t_h = 1/cosmo.H(0).decompose() 
flux = m_bh / (t_h * a_bh * m_p)
print("Hubble constant = {0:.2e}".format(cosmo.H(0)))
print("Hubble time = {0:.2e}".format(t_h))
print("flux = {0:.2e}".format(flux))

But dark energy (DE) can result from negatons through $P=-nmc^2$ (negative $\rho$ implying repulsive gravity). So density of such particles giving rise to DE can be estimated as follows. Assuming a separation between the particles of a Compton length, i.e. $\hbar/mc$, which is also positive for negatons (since both ℏ and m are negative) (negatons also obey the quantum uncertainty principle), the density is given by, $\frac{m}{(\hbar/mc)^3} =\frac{m^4 c^3}{\hbar ^3}$ , which is negative (since $m^4$ is positive and $\hbar^3$ is negative for negatons). Equating this to the observed DE density, $\rho_{DE}\approx 7\times 10^{-27} kgm^{-3}$, we can estimate the mass of these particles as:

$m=\left(\frac{\rho_{DE} \hbar^3}{c^3}\right)^{1/4}\approx 4\times 10^{-39} kg \quad (10)$

For these negative mass particles to make up dark energy of density $\approx 7\times 10^{-27} kgm^{-3}$, the number density is given by, $n \times m \approx 7\times 10^{-27} kgm^{-3}$. This gives:

$n\approx 10^{12} m^{-3} \quad (11)$

At earlier epochs (high $z$), the matter (and DM) densities would have been higher and as the Universe expands, the negative pressure due to the negatons starts to dominate, hence accounting for the accelerated expansion of the Universe observed at the present epoch.


rho_de = 7.e-27*u.kg*u.m**-3
m = ((rho_de * hbar**3) / c**3)**(1/4)
print("mass of negatons = {0:.2e}".format(m.decompose()))

n = rho_de / m
print("number density of negatons = {0:.2e}".format(n.decompose()))

##4. Consequences of negative mass 
Flux of negatons falling on a compact positive mass can reduce its mass. For instance a positive mass black hole (BH) accreting negatons could presumably end up having a negative mass. Here we have some consistency with cosmic censorship. 


### 4.1 Consequences for black holes
The effective mass of a Kerr-Newman BH is, $M_{eff}=M-a^2/r-Q^2/r$, so $M_{eff}$ can become negative if $a^2$, $Q^2$ reach extremal values. So also for charged BH. Extremal BH have $M^2=a^2+Q^2$, etc. By accreting negatons, the BH can become of effective negative mass. However once we have a negative mass BH it would not further capture any particles (both positive and negative) as a negative mass source repels all particles. Such a negative mass BH would also be described by the Schwarzschild solution (as both $G$,$M$ are reversed in sign). So it would have a horizon radius, $R=2GM/c^2$ . This would ensure symmetry of GR (and its solutions) under mass reversal, $M\rightarrow -M$. 

Can we distinguish a negative mass BH from a positive mass one? A negative mass BH cannot gather or accrete matter. So there is no quasar like activity. But it can still emit negative energy Hawking radiation. The Hawking temperature would be $T=(hc^3)/8\pi k_B GM$, but the quantum of action is negative for negatons. So it is seen $T$ is negative. Indeed to avoid a Boltzmann exponential blow up, negatons must have a negative temperature. However the entropy of a negative mass BH would be positive as, $S\propto (GM^2)/\hbar c$ ($G$,$\hbar$ negative, $M^2$ positive). This is consistent with second law of thermodynamics which cannot be violated. The evaporation lifetime would also be positive, given by $t=(G^2 M^3)/(\hbar c^4)$. With a negative T, black body emission energy is negative as $E=aT^4,a=k^5/(45\hbar^3 c^3)$, $k$ being positive and $\hbar$ negative. The entropy of the black body radiation is however positive ($\propto aT^3$, both negative). Again consistent with thermodynamics. 


# Calculate the Hawking temperature
T = lambda M: (h*c**3) / (8*np.pi*k_B*G*M)
n = np.linspace(24,40)
M = 10.**n*u.kg

# Calculate the evapouration timescale
t = lambda M: (G**2*M**3)/(hbar*c**4)


# Hawking tempertaure vs Mass
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5))
ax1.plot(n,np.log10(T(M).value),lw=2,c='b')
ax1.set_xlabel('log10 M (kg)')
ax1.set_ylabel('log10 T (K)')
ax1.grid(True)

# Evapouration timescale vs Mass
ax2.plot(n,np.log10(t(M).value),lw=2,c='r')
ax2.set_xlabel('log10 M (kg)')
ax2.set_ylabel('log10 t (s)')
ax2.grid(True)
plt.show()

### 4.2 Charge on negatons

If negatons have an electric charge, the charge squared, i.e. $e^2$ is negative. With this charge, negatons can form atoms (with the usual Bohr radii) but with negative energy levels. Also Larmor radiation energy would be negative (proportional to $a^2 e^2$, $a$ the acceleration) which implies that the energy would not be unbounded from below which is unphysical. 
In order to avoid reducing significantly the mass of SMBH or red supergiant stars (over a billion years) the limit background energy density in interstellar space of such negatons is easily shown to be $\rho <|10^{-25} g cm^{-3}|$, consistent with DM density. Like in an electrostatic plasma, a background density of both positons and negatons can give rise to a long range gravitational screening analogous to Debye screening. In Newtonian case, such screening can be described by a modified Poisson equation, 

$\nabla ^2 \phi - \frac{1}{\lambda_D^2}\phi=-4\pi G nm \delta(r) \quad (12)$

This will give rise to $\lambda_D \approx 2\pi \nu / (4 \pi Gnm)^{1/2}$ , which would be expected to be of the order of several kiloparsec, with $nm=\rho \approx 10^{-25} gcm^{-3}$. Such a modified potential (over kiloparsec distance) was used in [17, 18] to obtain flat rotation curves. Presence of such a distribution of negatons over megaparsec scales could prevent formation of galaxies, possibly accounting for presence of large voids free of galactic systems. 


## 5. Negaton bound states

Even negative mass can form bound states similar to Bohr's atom. The radius of the orbit is given as 

$r = \frac{4\pi\epsilon_0 n^2 \hbar^2}{m e^2} \quad (13)$.

Here the mass $m$  is negative and $\hbar^2$ is positive. So, for the radius to be positve, equation (13) implies $e^2$ is negative. Therefore the chrage on the negaton $e_{neg} = i e$. 

For the bound state we have the centripetal force balancing the Coloumb force ($F_{cen} = F_{Col}$) , i.e.,

$ \frac{m v^2}{r} = \frac{e^2}{4\pi \epsilon_0 r^2} \quad (14).$

For negatons, the negative sign of $m$ will cancel out the negative sign of $e^2$, hence we have the usual expression for kinetic energy from equation (14) as

$\frac{1}{2}mv^2= \frac{e^2}{8 \pi \epsilon_0 r} \quad (15)$.

The energy of the bound state is then given by 

$E = \frac{m e^4}{8\epsilon_0^2n^2h^2} \quad (16)$.

Here, $e^4$ is positive, $h^2$ is positive and $m$ is negative.

##6. Conclusion and future directions

Many recent papers have invoked negative mass as a possible model for dark energy. Here we have brought out some of the inconsistencies in these earlier work while considering the various gravitating properties of negative mass particles, and show that for a negative mass source (to have a positive binding energy associated with repulsive interaction) $G$ must change to $-G$. These discussions are carried over to the field equations in general relativity, making the field equations – for matter consisting of both positive and negative masses – consistent. We have also shown that the quantum vacuum of negatons, in the initial epoch, also gives rise to an accelerated expansion (inflationary scenario).

In short, the presence of negative masses (not ruled out by gravitational theories including GR) and their interesting gravitational interactions opens up a whole host of exciting possibilities for astrophysics, cosmology and expands our basic concepts in theoretical physics [17-19]. In particular, the implications for formation of black holes and their associated thermodynamics (including Hawking radiation), consequences for other compact objects (like white dwarfs and neutron stars), as well as for cosmic structure formation would be of much interest. 


## References:
[1] Jammer, M.: Concepts of Mass in Classical and Modern Physics. Harvard Univ. Press, Cambridge (1961)

[2] Bondi, H.: Negative mass in general relativity. Rev. Mod. Phys. 29, 423 (1957)

[3] Terletsky, Ya.P.: Positive, negative and imaginary rest masses. J. Phys. Radium, 23, 910 (1963)

[4] Farnes, J.S.: A unifying theory of dark energy and dark matter: Negative masses and matter creation within a modified ΛCDM framework. Astron. Astrophys. 620, A92 (2018)

[5] Petit, J.P., d’Agostini, G.: Negative mass hypothesis in cosmology and the nature of dark energy. Astrophys. Space Sci. 354, 611 (2014)

[6] Benoit-Lévy, A., Chardin, G.: Introducing the Dirac-Milne universe. Astron. Astrophys. 537, A78 (2012)

[7] Einstein, A.: Comment on Schrödinger’s Note “On a System of Solutions for the Generally Covariant Gravitational Field Equations”. Physikalische Zeitschrift, 19, 165 (1918). Translated by Engel A. in The Collected Papers of Albert Einstein, Vol. 7, The Berlin Years: Writings, 1918–1921 (Princeton University Press).

[8] Messiah, A.: Quantum mechanics. Interscience Publishers, New York (1961)

[9] Sakurai, J.J.: Mass reversal and weak interactions. Il Nuovo Cimento. 7, 649 (1958)

[10] Ashcroft, N.W., Mermin, N.D.: Solid State Physics. Harcourt, New York (1976)

[11] Morris, M.S. et al.: Wormholes, time machines, and the weak energy condition. Phys. Rev. Lett. 61, 1446 (1988)

[12] Hawking, S.W.: Particle creation by black holes. Comm. Math. Phys. 43, 199 (1975)

[13] Caldwell, R.R.: A phantom menace? Cosmological consequences of a dark energy component with super-negative equation of state. Phys. Lett. B. 545, 23 (2002) 

[14] Winterberg, F.: Remark concerning the gravitational interaction of matter and anti-matter. Il Nuovo Cimento. 19, 186 (1961)

[15] Forward, R.L.: Negative matter propulsion. J. Propul. Power. 6, 28 (1990)

[16] Buchdahl, H.A.: General relativistic fluid spheres. Phys. Rev. 116, 1027 (1959)

[17] de Sabbata, V., Sivaram, C.: Spin and Torsion in Gravitation. World Scientific, Singapore, (1994)

[18] Sivaram, C., Arun, K.: Some enigmatic aspects of the early universe. Astrophys. Space Sci. 334, 225 (2011)

[19] Sivaram, C., Arun, K.: Primordial rotation of the universe, hydrodynamics, vortices and angular momenta of celestial objects. Open Astron. J. 5, 7 (2012)
