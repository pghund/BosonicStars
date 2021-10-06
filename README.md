# BosonicStars

This is all code to try to compute the structure of a star with at least one bosonic species, the system of equations can be derived from Michael Kiessling's 2009 paper 
"Monotonicity of Quantum Ground State Energies: Bosonic Atoms and Stars". 

Boson2SpeciesEuler is just a standard Euler method ODE solver for a spherically symmetric star composed of only electrons and alpha particles. It does give a result that looks reasonable. But by a result of Lieb and Simon in the paper "The Hartree-Fock theory for Coulomb systems", any boson density must have exponential decay, so the results here cannot be correct. Why it does not work is unclear. 

BosonicSCF attempts to use the iterative solver of Hachisu, a description can be found in the repository on Fermionic stars, to find the densities. BosonicGradientDescent attempts to use gradient descent. A mathematical description of this approach can be found in the PDF BosonGradientDescent. 
