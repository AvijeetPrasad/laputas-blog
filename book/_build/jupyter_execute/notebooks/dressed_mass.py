<a href="https://colab.research.google.com/github/AvijeetPrasad/laputas/blob/main/notebooks/dressed_mass.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Dressed mass and singularity

## Summary:

In classical theories, as the size of the particle tends to zero, the energies associated with the particle approached infinity. Here we look at how interaction can "dress" a particle giving a total mass different than its bare mass. As a consequence, we see that the size of the particle becoming zero does not cause its mass to become infinite. The Mach's principle takes into account interactions of all the particles in the universe to give the inertial mass.
Using this idea, we can estimate the values of the interaction cross-sections of strong, weak, and electromagnetic forces. 

# Import relevant packages and constants
import math
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.constants.si import c, G, m_p, m_e, e, eps0, h, hbar
from astropy.cosmology import WMAP9 as cosmo
# c   = speed of light
# e   = charge of proton
# eps0= Vacuum electric permittivity
# G   =  universal gravatational constant
# h   = Planck's constant
# hbar= reduced Planck's constant
# m_e =  rest mass of electron
# m_p =  rest mass of proton

## 1. Introduction

The bare mass ($m_0$) of an elementary particle is the limit of its mass as the scale of distance approaches zero or, equivalently, as the energy of a particle collision approaches infinity. In general we can write the experimentally observed mass of a particles as: 

$m = m_0 + \delta_m \quad (1)$,

where $\delta_m$ is the change in mass due to interaction. 



## 2. Singularity

As the size of a particle tends to zero the energy associated with the particle tends to infinity. For a charge $e$ of size $r$, its electromagnetic self energy is  $e^2 / 4\pi \epsilon_0 r$. This tends to infinity as $r \rightarrow 0 $. Similarly for a particle of mass $m$ and size $r$, the gravitational self energy is given by $G m^2 / r$, which againtends to infinity as $r \rightarrow 0 $.

- For a particle of mass $m$ and charge $e$, the total mass can be written as:

  $m_{tot} = m_0 + \frac{e^2}{(4 \pi \epsilon_0  r c^2)} - \frac{G m_{tot}^2}{(r c^2)} \quad (2)$.

  Here $G m^2$ behaves as the gravitational charge. Solving for $m_{tot}$ as $r \rightarrow 0 $

  $m_{tot} = \frac{e}{\sqrt{4 \pi \epsilon_0 G}} \quad (3)$.

  This holds true even in General Relativity (ADM mass). 

- Considering rotational energy equation (2) becomes:

  $m_{tot} = m_0 + \frac{e^2}{(4\pi \epsilon_0 r c^2)} - \frac{G m_{tot}^2}{(r c^2)} + \frac{\hbar}{(rc)} \quad (4)$

  Solving for $m_{tot}$ as $r \rightarrow 0 $

  $m_{tot} = \sqrt{\frac{\hbar c}{G} + \frac{e^2}{4 \pi \epsilon_0 G}} \quad (5)$.

  From equations (3) and (5), we see that the bare mass does not appear in the expression as $r \rightarrow 0$.

  In classical model, as $r \rightarrow 0 $, the mass tends to infinity. Here, by including the self-energy terms, we have avoided this singularity.

- When gravity is as strong as the strong interaction (for which the coupling $\left(\frac{G m^2}{\hbar c} \sim 1 \right)$, we get:

  $m = \sqrt{\frac{\hbar c}{G}} \quad (6)$,

  which is same as the Plank mass ($m_{pl}$).

- With strong force $\left(\frac{g_s}{\hbar c} \approx 1 \right)$, a similar analysis gives:

  $m = \frac{g_s}{\sqrt{G}} \quad (7)$,

  where $g_s$ is the strong interaction charge.








# mtot1 = total mass without rotation
mtot1 = e / np.sqrt(4*np.pi*eps0*G)
mtot1 = mtot1.to(unit=u.kg)
print(f"Total mass without rotation = {mtot1:.2e}")

# mtot2 = total mass with rotation
mtot2 = np.sqrt((hbar * c / G) + (e**2/(4*np.pi*eps0*G)))
mtot2 = mtot2.to(unit=u.kg)
print(f"Total mass with rotation = {mtot2:.2e}")

### 2.1 Massive objects

Similar arguments can be applied for massive objects. If they have a charge $e$ or rotation $\hbar$, the mass will not tend to infinity as $r$ tends to zero. For example, in the case of the sun, with angular momentum $J$, as $r \rightarrow 0$, the mass is given by (similar to equation (6)):

$m = \sqrt{\frac{Jc}{G}} \quad (8)$.

This can be extended to all the matter in the universe collapsing to a zero size. The expression will remain the same even with relativistic mass ($\gamma m_0$).

### 2.2 Dressed mass

The "dressed" mass of a particle is the total mass including that due to the self energy.

- Particles with zero mass will have no charge, like photons for example, since they are not dressed by the interaction.

- The smallest mass with the fundamental charge is the electron with a mass of $0.511 MeV$.

- The proton which takes part in the strong interaction will be dressed by the strong force. Hence its mass is $\approx 1900$ times the electron mass.

- Since the weak interaction is about a million times weaker than elctromagentic the neutrino mass (dressed by the weak force) will be a million times smaller than the electron mass ($\sim 1 eV$).

## 3. Mach's principle

Mach’s principle says that a particle’s inertia is due to some interaction of that particle with all the other masses in the universe. The mass of a single isolated particle in the universe cannot be measured, we need at least two particles. According to Mach's principle all the inertial mass is due to the background interaction.

If the local mass $(m_{loc})$ is due to all the other masses ($M$) in the universe, and $r$ is the average radius of distant mass, then we have the equation:

$m_{loc} = \frac{G M m_{loc}}{r c^2} \quad (9)$,

where $M=\frac{4}{3}\pi\rho_{avg} r^3$, where $\rho_{avg}$ is the average density of the mass. Equation (9) implies $\frac{G M}{r c^2}\approx 1$. This could possibly explain why the value of $G$ is small, since $M$ is large.

If the average size $r = R_H$, the Hubble radius, the equation (9) now becomes

$8\pi G \rho_{av} T_H^2 \approx 1 \quad (10)$,

where $T_H$ is the Hubble time.


t_h = 1/cosmo.H(0).decompose() 
rho_av = 1/(8*math.pi*G*t_h**2)
print(f"The average density = {rho_av:0.2}")

The gravitational charge is given by $Gm^2$. The basic masses are the proton and the electron, hence the three possibilites are: $G m_p^2$, $G m_e^2$ and $G m_p m_e$. 
If $N$ is the total number of particles in the universe, then the fluctuations in these masses can be given by:

- $\sqrt{N} G m_p^2/\hbar c \sim 11.8 \quad (11a)$,

 which is analogous to the Pion constant (~14).

- $\sqrt{N} G m_e^2/\hbar c \sim 3.5 \times 10^{-6} \quad (11b)$,

 which is analogous to the Weak interaction constant.

- $\sqrt{N} G m_p m_e/\hbar c \sim 1/137 \quad (11c)$,

  which is analogous to the Fine-structure constant.

# N: The total number of particles in the Universe
N = 4.e78 
# fpi : Pion constant ~ 14.
fpi = ((np.sqrt(N) * G * m_p**2) / (hbar * c)).decompose() 
# fweak: Weak interaction constant is between 10^-6 and 10^-7
fweak = ((np.sqrt(N) * G * m_e**2) / (hbar * c)).decompose() 
# alpinv : Fine structure constant = 1 / 137 
alpinv = ((np.sqrt(N) * G * m_e*m_p) / (hbar * c)).decompose()  

print(f"The interaction cross-sections are: \n\
strong: {fpi:0.2}\nweak  : {fweak:0.2e}\nEM    :{alpinv:0.2}")