# Serpentine
This Matlab code describes the electrical and mechanical behaviors of the SCARS sensor, with input variables of the beam geometry, prestretch value, and the desired resolution for strain-resistance calculation

mech2res_num is a numerical solution that takes the serpentine geometry and amount of prestretch to find the beam position that results in minimum energy stored. The energy is composed of two components: the strain energy in the beam, which is found through numerical integration, and the strain energy in the substrate, which is found by integrating along each beam length given the beam's displacement from its equilibrium (i.e., zero substrate energy) position. 

mech2res_closed uses a closed form solution to find the displacement of each beam by first assuming the displacement of the beam's tip is either the displacement to its next neighbor (displacement of gap -a) or the displacement that minimizes the energy stored in the substrate for all following (right neighbor +) beams. The function then finds the "effective length" of the beam that causes the minimum energy displacement at the end of the effective length. The remaining length of the beam beyond the effective length is assumed to mirror the previous (left neighbor) beam and has no stored beam strain energy.

Within the mech2res_closed function, the code may be changed to solve for a constant distributed force along the beam (dU_vec_dist) or a point force at the end of the beam (dU_vec_point). 

Each mech2res* uses the same area_to_resistance function to find the total network resistance. This function takes in the contact length of each beam and uses it to solve for the contact resistance between each beam as well as the total resistance of the serpentine. 

