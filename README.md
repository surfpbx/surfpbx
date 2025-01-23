# Surface pourbaix diagram generation utilities

This repository contains scripts for generating surface Pourbaix diagrams and acqueous stability slices. Examples on how to generate and display diagrams are provided as well.

If you use this code or the ASE standard Pourbaix implementation, please cite
"Predicting aqueous and electrochemical stability of 2D materials from extended Pourbaix analyses" - S. Americo, I.E. Castelli and K.S. Thygesen

The "data" directory contains the complete set of materials discussed in the publication.
The material set is presented as a text file containing the material formulas, their formation energies, their convex hull energies,
and Pourbaix energies in different electrochemical conditions.

Requires ASE >= 3.24.0.
