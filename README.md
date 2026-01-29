# Supernova quantities modelization

Project Status: Completed (Academic Project)
Development Period: Nov 2021 – Jan 2022
Language: C

## Description
This is a standalone C implementation of a numerical solver describing a strong spherical blast wave propagating in a uniform medium. The solver computes the dimensionless profiles of density, velocity, and pressure behind a shock front, corresponding to the Sedov–Taylor solution for an adiabatic explosion in a polytropic gas. It models the physics of Supernova using the Runge-Kutta 4th order method.

## How to use
The code is contained in a single file for ease of use.
Compile with any standard C compiler:
$ gcc rico.c -o rico.x -lm
