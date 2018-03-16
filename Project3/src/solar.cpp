#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>

using namespace std; 

const double pi=3.14159265359; 
const double Gconst=4*pi*pi; 

int main(int argc, char* argv[])
{
    // command line input
    string filename; 
    double time,dt; 
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

    //Sun-Earth system 
    double earth_r[3]={0},earth_v[3]={0},earth_a[3]={0}; //0, 1, 2 for x, y, z direction, respectively
    double sun_r[3]={0},sun_v[3]={0},sun_a[3]={0};
    double earth_force[3]={0},sun_force[3]={0}; 
    double earth_mass, sun_mass; 
    //initialize
    earth_r[0]=1.0; earth_r[1]=0.0; earth_r[2]=0.0; 
    earth_v[0]=0.0; earth_v[1]=2*pi; earth_v[3]=0.0; 
    earth_mass=3e-6; sun_mass=1.0; 
    ofstream outfile; 
    double distsqr; 
    //Euler's method
    outfile.open(filename+"_Euler.txt"); 
    outfile << "Earth's x y z; Sun's x y z" <<endl; 
    int n; 
    n=int(time/dt); 
 
    for (int i=0; i<=n; i++)
    {
        distsqr=0.0; 
        for (int j=0; j<3; j++)
            distsqr=distsqr+(earth_r[j]-sun_r[j])*(earth_r[j]-sun_r[j]);
        for (int j=0; j<3; j++)
        {
            earth_force[j]=Gconst*sun_mass*earth_mass/distsqr*(sun_r[j]-earth_r[j])/sqrt(distsqr); 
            earth_a[j]=earth_force[j]/earth_mass; 
            earth_v[j]=earth_v[j]+earth_a[j]*dt; 
            earth_r[j]=earth_r[j]+earth_v[j]*dt;  
            sun_force[j]=-earth_force[j]; 
            sun_a[j]=sun_force[j]/sun_mass; 
            sun_v[j]=sun_v[j]+sun_a[j]*dt; 
            sun_r[j]=sun_r[j]+sun_v[j]*dt;  
        }
        for (int j=0; j<3; j++)
            outfile << earth_r[j] <<' '; 
        for (int j=0; j<3 ;j++)
            outfile << sun_r[j] <<' '; 
        outfile <<endl; 
    }
    outfile.close(); 
    
    
    //initialize
    earth_r[0]=1.0; earth_r[1]=0.0; earth_r[2]=0.0; 
    earth_v[0]=0.0; earth_v[1]=2*pi; earth_v[3]=0.0; 
    sun_r[0]=0.0; sun_r[1]=0.0; sun_r[2]=0.0; 
    sun_v[0]=0.0; sun_v[1]=0.0; sun_v[2]=0.0; 
    earth_mass=3e-6; sun_mass=1.0; 
    double earth_olda[3]={0},sun_olda[3]={0}; 
    //Verlet method
    outfile.open(filename+"_Verlet.txt"); 
    outfile << "Earth's x y z; Sun's x y z" <<endl; 
    n=int(time/dt); 
    distsqr=0.0; 
    for (int j=0; j<3; j++)
        distsqr=distsqr+(earth_r[j]-sun_r[j])*(earth_r[j]-sun_r[j]); 
    for (int j=0; j<3; j++)
    {
        earth_force[j]=Gconst*sun_mass*earth_mass/distsqr*(sun_r[j]-earth_r[j])/sqrt(distsqr); 
        sun_force[j]=-earth_force[j]; 
        earth_a[j]=earth_force[j]/earth_mass; 
        sun_a[j]=sun_force[j]/sun_mass;
    }
    for (int i=0; i<=n; i++)
    {
        for (int j=0; j<3; j++)
        {            
            earth_r[j]=earth_r[j]+dt*earth_v[j]+0.5*dt*dt*earth_a[j]; 
            sun_force[j]=-earth_force[j]; 
            sun_r[j]=sun_r[j]+dt*sun_v[j]+0.5*dt*dt*sun_a[j]; 
            earth_olda[j]=earth_a[j]; 
            sun_olda[j]=sun_a[j]; 
        }
        
        distsqr=0.0; 
        for (int j=0; j<3; j++)
            distsqr=distsqr+(earth_r[j]-sun_r[j])*(earth_r[j]-sun_r[j]); 
        for (int j=0; j<3; j++)
        {
            earth_force[j]=Gconst*sun_mass*earth_mass/distsqr*(sun_r[j]-earth_r[j])/sqrt(distsqr); 
            sun_force[j]=-earth_force[j]; 
            earth_a[j]=earth_force[j]/earth_mass; 
            sun_a[j]=sun_force[j]/sun_mass; 
            earth_v[j]=earth_v[j]+0.5*dt*(earth_a[j]+earth_olda[j]); 
            sun_v[j]=sun_v[j]+0.5*dt*(sun_a[j]+sun_olda[j]);             
        }
        
        for (int j=0; j<3; j++)
            outfile << earth_r[j] <<' '; 
        for (int j=0; j<3 ;j++)
            outfile << sun_r[j] <<' '; 
        outfile <<endl; 
    }
    outfile.close(); 
    
    return 0; 
}