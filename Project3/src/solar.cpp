#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "ode.hpp"
#include "planet.hpp"

using namespace std; 

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
    if (!infile) 
    {
        cerr << "Cannot open input file!"; 
        return 1; 
    }
    int num; 
    infile >>method>>num; 
    planet **solar=NULL; 
    if (num<=0) return 0; 
    solar=new planet*[num]; 
    if (solar==NULL) 
    {
        cerr << "Cannot allocate memory for 'solar'!"; 
        return 1; 
    }
    string readname; 
    double ma,x0,y0,z0,vx0,vy0,vz0; 
    bool fix,comcorrect; 
    for (int i=0; i<num; i++)
    {
        solar[i]=NULL; 
        infile >>readname>>ma>>x0>>y0>>z0>>vx0>>vy0>>vz0>>fix; 
        solar[i]=new planet(readname,ma,x0,y0,z0,vx0,vy0,vz0,fix); 
        if (solar[i]==NULL)
        {
            cerr << "Cannot allocate memory for a planet!"; 
            return 1; 
        }
    }
    infile >>comcorrect; //whether move the center of mass to origin 
    if (comcorrect) com_correct(solar, num); 
    
    int n; 
    n=int(time/dt); 
    
    ofstream* outfile=NULL; 
    outfile=new ofstream[num]; 
    if (outfile==NULL)
    {
        cerr << "Cannot allocate memory for 'outfile' array"; 
        return 1; 
    }
    for (int i=0; i<num; i++)
    {
        outfile[i].open(filename+"_"+to_string(i)+".txt"); 
        if (!outfile) 
        {
            cerr << "Cannot open output file!"; 
            return 1; 
        }
        outfile[i] <<"Time: "<<time<<endl<<"Time step: "<<dt<<endl<<"Method: "<<method<<endl;
        outfile[i] <<"time name mass x y z vx vy vz Ek Ep E lx ly lz Unbound?"<<endl;  
    }
    
    //start calculation
    if (method) 
        verlet(solar,num,dt,n,outfile);  //Verlet method 
    else
        euler(solar,num,dt,n,outfile);   //Euler's method
    
    for (int i=0; i<num; i++)
    {
        delete solar[i]; 
        outfile[i].close(); 
    }
    
    delete[] solar; 
    delete[] outfile; 
    
    return 0; 
}