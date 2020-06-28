# Blast-wave invariant yield and v2
Particle transvers momentum (pT) spectra and azimuthal anisotropy in ultra-relativistic heavy-ion collisions can be rather well described by a superposition of moving sources where in the rest frame of each source the particles follow a thermal Boltzmann-Jüttner distribution.

In [arXiv:1910.14618](https://arxiv.org/abs/1910.14618) (Klaus Reygers, Alexander Schmah, Anastasia Berdnikova, Xu Sun, [Phys. Rev. C 101, 064905](https://doi.org/10.1103/PhysRevC.101.064905) (2020)) compact formulas were derived which consistently describe pT spectra and the elliptic flow v2(pT).

This code implements Eq.(4) (invariant yield) and Eq.(9) (elliptic flow v2(pT)) of the above paper.    

The formulas are based on an elliptical freeze-out hyper-surface. It is assumed that the freeze-out time of a fluid cell is independent of the radial coordinate. The formula (Eq. 4) for the invariant yield contains the classic formula obtained by Schnedermann, Sollfrank, and Heinz [nucl-th/9307020](https://arxiv.org/abs/nucl-th/9307020), ([Phys. Rev. C 48, 2462 (1993)](https://doi.org/10.1103/PhysRevC.48.2462)) as a special case (spherical freeze-out surface corrosponding to head-on collisions of spherical nuclei)

The code uses the root framework (https://root.cern.ch/)

Run it with

	root blastwave_yield_and_v2_sudden_freeze_out_in_r.C
