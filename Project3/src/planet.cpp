#include <fstream>
#include <string>
#include <cmath>
#include "planet.hpp"

using namespace std; 

planet:: planet()
{
    name="unknown"; 
}

planet:: planet(const string & nn)
{
    name=nn; 
}

planet:: planet (const string & nn, double m)
{
    name=nn; 
    mass=m; 
}

planet:: planet (const string & nn, double m, double pos[3], double vel[3])
{
    name=nn; mass=m; 
    r[0]=pos[0]; r[1]=pos[1]; r[2]=pos[2]; 
    v[0]=vel[0]; v[1]=vel[1]; v[2]=vel[2]; 
    distance=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
}

planet:: planet (const string & nn, double m, double x, double y, double z, double vx, double vy, double vz)
{
    name=nn; mass=m; 
    r[0]=x; r[1]=y; r[2]=z; 
    v[0]=vx; v[1]=vy; v[2]=vz; 
    distance=sqrt(x*x+y*y+z*z); 
}

planet:: planet (const string & nn, double m, double pos[3], double vel[3],bool fix)
{
    name=nn; mass=m; 
    r[0]=pos[0]; r[1]=pos[1]; r[2]=pos[2]; 
    v[0]=vel[0]; v[1]=vel[1]; v[2]=vel[2]; 
    fixed=fix; 
    distance=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
}

planet:: planet (const string & nn, double m, double x, double y, double z, double vx, double vy, double vz, bool fix)
{
    name=nn; mass=m; 
    r[0]=x; r[1]=y; r[2]=z; 
    v[0]=vx; v[1]=vy; v[2]=vz; 
    fixed=fix; 
    distance=sqrt(x*x+y*y+z*z); 
}

double planet:: dist(planet & partner) const
{
    double d;
    d=0.0;     
    for (int i=0;i<3;i++)
        d=d+(r[i]-partner.r[i])*(r[i]-partner.r[i]);  
    return sqrt(d); 
}

void planet:: calc_force(planet & partner, double &fx, double &fy, double &fz) const
{
    double temp,di; 
    di=dist(partner); 
    temp=Gconst*mass*partner.m()/di/di/di; 
    fx=temp*(partner.r[0]-r[0]); 
    fy=temp*(partner.r[1]-r[1]); 
    fz=temp*(partner.r[2]-r[2]); 
}

void planet:: init_newstep()
{
    force[0]=force[1]=force[2]=0.0; 
    init=true; 
}

void planet:: force_add_both(planet & partner)
{
    double fx,fy,fz; 
    if (!init) init_newstep(); 
    if (!(partner.init)) partner.init_newstep(); 
    calc_force(partner,fx,fy,fz); 
    force[0]+=fx; force[1]+=fy; force[2]+=fz; 
    for (int i=0;i<3;i++)
        partner.force[i]+=-force[i]; 
}

void planet:: accelerate_update()
{
    if (!fixed) 
        for (int i=0; i<3; i++)
        {
            olda[i]=a[i]; 
            a[i]=force[i]/mass; 
        }
    init=false; 
}

void planet:: Euler_update(double dt)
{
    accelerate_update(); 
    if (!fixed)
    {
        for (int j=0; j<3; j++)
        {
            v[j]=v[j]+a[j]*dt; 
            r[j]=r[j]+v[j]*dt;  
        }
        distance_update(); 
    }
    time=time+dt; 
}

void planet:: Verlet_r(double dt)
{
    if (fixed) return; 
    for (int j=0; j<3; j++)
        r[j]=r[j]+dt*v[j]+0.5*dt*dt*a[j]; 
}

void planet:: Verlet_v(double dt)
{   
    accelerate_update(); 
    if (!fixed)
    {           
        for (int j=0; j<3; j++)
            v[j]=v[j]+0.5*dt*(a[j]+olda[j]); 
        distance_update(); 
    }
    time=time+dt; 
}


void planet:: fileoutput(ofstream &file) const
{
    file <<name<<' '<<mass<<' '<<r[0]<<' '<<r[1]<<' '<<r[2]<<' '<<v[0]<<' '<<v[1]<<' '<<v[2]<<' '<<monofar<<endl; 
}

void planet:: distance_update()
{
    double temp; 
    temp=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    if ((monofar) && (temp<distance)) monofar=false; 
    distance=temp; 
}