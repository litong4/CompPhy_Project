#ifndef ode_H
#define ode_H

#include <fstream>
#include "planet.hpp" 

void calc_all_force(planet **, int);
void calc_all_potential(planet **,int); 
void verlet(planet **, int , double, int, ofstream*); //Verlet method
void euler(planet **, int , double, int, ofstream*); //Euler's method

#endif
