# README

This is documentation detailing the research conducted for acceleration noise for the IIP grant awarded to PSSL at the University of Florida.

There is a compilation of files and data, some used and some not, but here are all of the things that I have worked on.

## Comprehensive Acceleration noise

The comprehensive plot file is titled paperaccelnoise.m. For the magnetic accleration noise, there is data that has come from [https://podaac.jpl.nasa.gov/](https://podaac.jpl.nasa.gov/) - these files are stored in GraceFO Data. To grab another day's data, refer to [https://podaac.jpl.nasa.gov/](https://podaac.jpl.nasa.gov/). 

## Molflow Model

The data from running molflow simulations is found in Molflow Simulations.
To run a molflow simulation, load up the geometry by loading the autosave file located in /Molflow Simulations/molflow_win_2.7.10/molflow_win_2.7.10. Use Outgassing Coefficient.xlsx calculator to calculate the outgassing coefficient for the desired temperature, and then modify the desorption of the Molflow geometry to suit the current one.