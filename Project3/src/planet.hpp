#ifndef planet_H
#define planet_H

#include <iostream>
#include <fstream> 
#include <string>
#include <cmath>

using namespace std; 

const double pi=3.14159265359; 
const double Gconst=4*pi*pi; 

class planet
{
private:
    double r[3]={0},v[3]={0},a[3]={0},force[3]={0},olda[3]={0}; //0, 1, 2 for x, y, z direction, respectively
    double mass=0.0; 
    string name; 
    double time=0.0;
    bool init=false; 
    double dist(planet & partner) const; 
    void calc_force(planet & partner, double &fx, double &fy, double &fz) const; 
    void init_newstep(); 
public: 
    planet (); 
    planet (const string & nn); 
    planet (const string & nn, double m); 
    planet (const string & nn, double m, double pos[3], double vel[3]); 
    planet (const string & nn, double m, double x, double y, double z, double vx, double vy, double vz); 
    double t() const {return time;}
    double ax() const {return a[0];}
    double ay() const {return a[1];}
    double az() const {return a[2];}
    double vx() const {return v[0];}
    double vy() const {return v[1];}
    double yz() const {return v[2];}
    double rx() const {return r[0];}
    double ry() const {return r[1];}
    double rz() const {return r[2];}
    double rr(int i) const {if ((i<0)||(i>2)) return 0; else return r[i];}
    double vv(int i) const {if ((i<0)||(i>2)) return 0; else return v[i];} 
    double aa(int i) const {if ((i<0)||(i>2)) return 0; else return a[i];}
    string callname() const {return name;}
    double m() const {return mass;}
    void force_add_both(planet & partner); 
    void accelerate_update(); 
    void Euler_update(double dt); 
    void Verlet_r(double dt); 
    void Verlet_v(double dt); 
    void fileoutput(ofstream &file) const; 
};


#endif
