# Serpentine
Matlab code to describe electrical and mechanical behavior of serpentine

mech2resadjstrain is a numerical solution that takes the serpentine geometry and amount of prestretch to find the beam position that results in minimum energy stored. The energy is composed of two components: the strain energy in the beam, which is found through numerical integration, and the strain energy in the substrate, which is found by integrating along each beam length given the beam's displacement from its equilibrium (i.e., zero substrate energy) position. 

Each mech2res* uses the same area to resistance length function to find the total network resistance. This function takes in the overlap length of each beam and uses it to solve for the contact resistance between each beam as well as the total resistance of the serpentine. 

Area_to_res_length function takes the effective length of each beam into account when calculating "Rs" and "Rl"
