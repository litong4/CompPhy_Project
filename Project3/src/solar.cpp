#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "planet.hpp"

using namespace std; 

inline void calc_all_force(planet **solar,int num)
{
    for (int j=0; j<num; j++)
        for (int k=j+1; k<num; k++)
            solar[j]->force_add_both(*(solar[k])); 
}

int main(int argc, char* argv[])
{
    // command line input
    string filename; 
    double time,dt; 
    int method; //0: Euler. other: Verlet
    if (argc != 4)
    {
        cout << "Error in the input arguments! Please input two arguments: time t, time step dt and output filename. "<<endl;
        return 0; 
    }
    else
    {
        time=atof(argv[1]);
        dt=atof(argv[2]); 
        filename=argv[3];
    }

    //solar system initialize
    ifstream infile; 
    infile.open(filename+"_input.txt"); 
    int num; 
    infile >>method>>num; 
    planet **solar; 
    if (num<=0) return 0; 
    solar=new planet*[num]; 
    string readname; 
    double ma,x0,y0,z0,vx0,vy0,vz0; 
    for (int i=0; i<num; i++)
    {
        infile >>readname>>ma>>x0>>y0>>z0>>vx0>>vy0>>vz0; 
        solar[i]=new planet(readname,ma,x0,y0,z0,vx0,vy0,vz0); 
    }
    
    ofstream outfile; 
    outfile.open(filename+"_output.txt"); 
    int n; 
    n=int(time/dt); 
    outfile <<"Time: "<<time<<endl<<"Time step: "<<dt<<endl<<"Method: "<<method<<endl; 
    //start calculation
    if (method) 
    {//Verlet method
        calc_all_force(solar,num); 
        for (int j=0; j<num; j++) solar[j]->accelerate_update(); 
        for (int i=0; i<=n; i++)
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
    else
    {//Euler's method
        for (int i=0; i<=n; i++)
        {
            calc_all_force(solar,num); 
            for (int j=0; j<num; j++)
            {
                solar[j]->Euler_update(dt); 
                solar[j]->fileoutput(outfile); 
            }
        }
    }
    outfile.close(); 
    
    return 0; 
}