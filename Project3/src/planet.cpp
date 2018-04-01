#include <iostream>
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
    kinetic=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    ang_update(); 
}

planet:: planet (const string & nn, double m, double x, double y, double z, double vx, double vy, double vz)
{
    name=nn; mass=m; 
    r[0]=x; r[1]=y; r[2]=z; 
    v[0]=vx; v[1]=vy; v[2]=vz; 
    distance=sqrt(x*x+y*y+z*z); 
    kinetic=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    ang_update(); 
}

planet:: planet (const string & nn, double m, double pos[3], double vel[3],bool fix)
{
    name=nn; mass=m; 
    r[0]=pos[0]; r[1]=pos[1]; r[2]=pos[2]; 
    v[0]=vel[0]; v[1]=vel[1]; v[2]=vel[2]; 
    fixed=fix; 
    distance=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    kinetic=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    if ((fix)&&(kinetic>1e-10)) 
    {
        cerr << "Fixed planet should not have any kinetic energy!"; 
        exit(2); 
    }
    ang_update(); 
}

planet:: planet (const string & nn, double m, double x, double y, double z, double vx, double vy, double vz, bool fix)
{
    name=nn; mass=m; 
    r[0]=x; r[1]=y; r[2]=z; 
    v[0]=vx; v[1]=vy; v[2]=vz; 
    fixed=fix; 
    distance=sqrt(x*x+y*y+z*z); 
    kinetic=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    if ((fix)&&(kinetic>1e-10)) 
    {
        cerr << "Fixed planet should not have any kinetic energy!"; 
        exit(2); 
    }
    ang_update(); 
}

double planet:: dist(planet & partner) const
{
    double d;
    d=0.0;     
    for (int i=0;i<3;i++)
        d=d+(r[i]-partner.r[i])*(r[i]-partner.r[i]);  
    return sqrt(d); 
}

void planet:: calc_force(planet & partner, double &fx, double &fy, double &fz, double &pot) const
{
    double temp,di; 
    di=dist(partner); 
    temp=Gconst*mass*partner.m()/di/di/di; 
    fx=temp*(partner.r[0]-r[0]); 
    fy=temp*(partner.r[1]-r[1]); 
    fz=temp*(partner.r[2]-r[2]); 
}

double planet:: calc_pot(planet &partner) const
{
    return (-Gconst*mass*partner.m()/dist(partner)); 
}


void planet:: init_newstep()
{
    force[0]=force[1]=force[2]=0.0; 
    init=true; 
}

void planet:: init_potcal()
{
    potential=0.0; 
    init_pot=true; 
}

void planet:: force_add_both(planet & partner)
{
    double fx,fy,fz,pot; 
    if (!init) init_newstep(); 
    if (!(partner.init)) partner.init_newstep(); 
    calc_force(partner,fx,fy,fz,pot); 
    force[0]+=fx; force[1]+=fy; force[2]+=fz; 
    for (int i=0;i<3;i++)
        partner.force[i]+=-force[i]; 
}

void planet:: pot_add_both(planet & partner)
{
    double pot; 
    if (!init_pot) init_potcal(); 
    if (!(partner.init_pot)) partner.init_potcal(); 
    pot=calc_pot(partner); 
    potential+=pot; 
    partner.potential+=pot; 
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

void planet:: ang_update()
{
    for (int j=0; j<3; j++)
        ang[j]=mass*(r[(j+1)%3]*v[(j+2)%3]-r[(j+2)%3]*v[(j+1)%3]); 
}

void planet:: distance_update()
{
    double temp; 
    temp=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    if ((monofar) && (temp<distance)) monofar=false; 
    distance=temp; 
}

void planet:: Euler_update(double dt)
{
    accelerate_update(); 
    kinetic=0.0; 
    if (!fixed)
    {
        for (int j=0; j<3; j++)
        {
            v[j]=v[j]+a[j]*dt; 
            r[j]=r[j]+v[j]*dt;  
            kinetic+=0.5*mass*v[j]*v[j]; 
        }
        ang_update(); 
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
    kinetic=0.0; 
    if (!fixed)
    {           
        for (int j=0; j<3; j++)
        {
            v[j]=v[j]+0.5*dt*(a[j]+olda[j]); 
            kinetic+=0.5*mass*v[j]*v[j]; 
        }
        ang_update(); 
        distance_update(); 
    }
    time=time+dt; 
}


void planet:: fileoutput(ofstream &file) const
{
    file <<time<<' '<<name<<' '<<mass<<' '
    <<r[0]<<' '<<r[1]<<' '<<r[2]<<' '
    <<v[0]<<' '<<v[1]<<' '<<v[2]<<' '
    <<kinetic<<' '<<potential<<' '<<energy<<' '
    <<ang[0]<<' '<<ang[1]<<' '<<ang[2]<<' '<<monofar<<endl; 
}

void planet:: energy_update()
{
    energy=kinetic+potential; 
    init_pot=false; 
}

void com_correct(planet **solar, int num) //change to center of mass frame 
{
    double totalmass; 
    double rcm[3]={0},vcm[3]={0}; 
    totalmass=0.0; 
    for (int i=0; i<num; i++)
    {
        totalmass+=solar[i]->m(); 
        for (int j=0; j<3; j++)
        {
            rcm[j]+=(solar[i]->m())*(solar[i]->rr(j)); 
            vcm[j]+=(solar[i]->m())*(solar[i]->vv(j)); 
        }
    }
    for (int j=0; j<3; j++)
    {
        rcm[j]=rcm[j]/totalmass; 
        vcm[j]=vcm[j]/totalmass; 
    }
    for (int i=0; i<num; i++)
        for (int j=0; j<3; j++)
        {
            solar[i]->r[j]=(solar[i]->r[j])-rcm[j]; 
            solar[i]->v[j]=(solar[i]->v[j])-vcm[j]; 
        }
    /* //check whether we successfully change to c.o.m frame 
    double checkp[3]={0}, checkr[3]={0}; 
    for (int i=0; i<num; i++)
        for (int j=0; j<3; j++)
        {
            checkp[j]+=(solar[i]->mass)*(solar[i]->v[j]); 
            checkr[j]+=(solar[i]->mass)*(solar[i]->r[j]); 
        }
    cout <<checkp[0]<<' '<<checkp[1]<<' '<<checkp[2]<<endl; 
    cout <<checkr[0]<<' '<<checkr[1]<<' '<<checkr[2]<<endl; 
    */
}
