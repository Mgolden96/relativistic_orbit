# relativistic_orbit
This code numerically integrates orbits around a Schwarzschild black hole. It is assumed throughout that G=c=M=1.

To use, simply compile main.cpp with your favorite C++ compiler and mess around with parameters.dat. Use your favorite plotting software to see the resulting orbits. The output of the simulation is dumped to orbit.dat.

This code uses third-order Runge-Kutta and Schwarzschild coordinates for the spacetime. Evolution is carried out with respect to proper time, not coordinate time.
