Read Me — CAI Condensation/Evaporation Model

This code accompanies the manuscript “Ca isotopes in CAIs indicate rapid condensation of solids in the early solar nebula” submitted to EPSL in 2025.  

The code is being shared for transparency and archiving purposes; it was written to perform specific calculations for a manuscript.  It was not designed to be used-friendly or robust for all purposes.  Time permitting, it will be updated and reorganized in the future to be cleaner and more useful / user-friendly.

Files/folders:

- CAI_EvapCond_MainScript_EPSL:
This is the main script to run the model.  It includes an “inputs” section at the top, where input values can be adjusted by the user.  [Note, however, that this code is being shared for transparency and archiving purposes; the code has not yet been designed to be entirely user-friendly, robust to all inputs, or valid over all conditions.]  At the very end of the script, there is a section that generates a plot.  

- functions/
Folder that contains multiple .m files that are functions relied upon by CAI_EvapCond_MainScript_EPSL.m

- data/
Folder that contains input thermodynamic values, elemental abundances, and other values relied upon by CAI_EvapCond_MainScript_EPSL

Instructions:

(1) Open and run the script (CAI_EvapCond_MainScript_EPSL) in MATLAB.  The folders “functions” and “data” should be located in the same folder as the main script.  

(2) Some important output variables that can be useful for plotting results include: 
	“t”, a vector with time values (in seconds); 
	“elements”, a list of elements, ordered as in below elemental arrays;
	“species”, a list of species, ordered as in below species array;
	“b_gss_array”, an array with abundances of each element in gas over all times;
	“b_css_array”, an array with abundances of each element in particle over all times;
	“b_isotopes_gss_array”, an array with abundances of isotopes in gas over all times;
	“b_isotopes_css_array”, an array with abundances of isotopes in particle over all times;
	“n_species_array”, an array with abundances of each species in gas and particle over all times.

(3) Separate plotting scripts are not included, but a simple example plot is generated for Al/Ca and d44Ca at the end of CAI_EvapCond_MainScript_EPSL.m.

Notes:

The code has been tested with MATLAB version 2023a.  

Among multiple limitations of the code: (A) The code and data files include only select gas and mineral species. So, if running at significantly non-solar compositions (e.g., different C/O or fO2), important species may be missing.  (B) There can sometimes be failures when running calculations with too small or large a step size; effort will likely be made in future to fix this.  (C) If melting occurs, evaporation of elements other than those in melt model (Ca, Mg, Si, Al), such as Ti, will not be properly calculated
