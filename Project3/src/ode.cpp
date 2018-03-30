#include <fstream>
#include "planet.hpp" 
#include "ode.hpp"

inline void calc_all_force(planet **solar,int num)
{
    for (int j=0; j<num; j++)
        for (int k=j+1; k<num; k++)
            solar[j]->force_add_both(*(solar[k])); 
}

void verlet(planet **solar, int num, double dt, int stepnum, ofstream &outfile) //Verlet method
{
    calc_all_force(solar,num); 
    for (int j=0; j<num; j++) solar[j]->accelerate_update(); 
    for (int i=0; i<stepnum; i++)
    {
        for (int j=0; j<num; j++) solar[j]->Verlet_r(dt); 
        calc_all_force(solar,num); 
        for (int j=0; j<num; j++) 
        {
            solar[j]->Verlet_v(dt); 
            solar[j]->fileoutput(outfile); 
        }
    }
}

void euler(planet **solar, int num, double dt, int stepnum, ofstream &outfile) //Euler's method
{
    //Euler's method
    for (int i=0;i<stepnum;i++)
    {
        calc_all_force(solar,num); 
        for (int j=0; j<num; j++)
        {
            solar[j]->Euler_update(dt); 
            solar[j]->fileoutput(outfile); 
        }
    }
}
