# Surface pourbaix diagram generation utilities

This repository contains the code used for generating surface Pourbaix diagrams and acqueous stability slices.

If you use this code or the ASE standard Pourbaix implementation, please cite
"Predicting aqueous and electrochemical stability of 2D materials from extended Pourbaix analyses" - S. Americo, I.E. Castelli and K.S. Thygesen

The scripts in the "examples" directory show how to use the code, reproducing Fig. 2b and [INSERT] in the publication.

The "data" directory contains the complete set of 3376 materials discussed in the conventional Pourbaix diagrams (CPD) section of the publication.
The material set is presented as a text file containing the material formulas, their convex hull energies, their formation energies, 
their Pourbaix energies in different electrochemical conditions.

**Code requirements:** ASE >= 3.24.0.
