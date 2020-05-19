# Hybrid-Beamforming

This example introduces the basic concept of hybrid beamforming and shows how to split the precoding and combining weights using orthogonal matching pursuit algorithm. It shows that hybrid beamforming can closely match the performance offered by optimal digital weights.
This program simulates a 64 x 16 MIMO hybrid beamforming system, with a 64-element square array with 4 RF chains on the transmitter side and a 16-element square array with 4 RF chains on the receiver side. Each antenna element can be connected to one or more TR modules.
It is assumed that each antenna is connected to all RF chains. Thus, each antenna is connected to 4 phase shifters. 
Such an array can be modeled by partitioning the array aperture into 4 completely connected subarrays.

Spectral Efficiency Comparison : 
Compares the spectral efficiency achieved using the optimal weights with that of the proposed hybrid beamforming weights. 
The resulting spectral efficiency curve is obtained from 50 Monte-Carlo trials for each SNR.
The transmit antenna array is assumed to be at a base station, with a focused beamwidth of 60 degrees in azimuth and 20 degrees in elevation. 
The hybrid beamforming can perform close to what optimal weights can offer using less hardware.
