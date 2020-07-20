# Serpentine
This Matlab code describes the electrical and mechanical behaviors of the SCARS sensor, with input variables of the beam geometry, prestretch value, and the desired resolution for strain-resistance calculation

mech2res_closed uses a closed form solution to find the displacement of each beam by first assuming the displacement of the beam's tip is either the displacement to its next neighbor (displacement of gap -a) or the displacement that minimizes the energy stored in the substrate for all following (right neighbor +) beams. The function then finds the "effective length" of the beam that causes the minimum energy displacement at the end of the effective length. The remaining length of the beam beyond the effective length is assumed to mirror the previous (left neighbor) beam. The resistance is calculated using this overlap length.

area_to_resistance is the electrical component of the model. It takes the contact area of each beam from its neighbors and solves for the thevenin resistance of the sensor using Nodal analysis with a test current of 1 A.

