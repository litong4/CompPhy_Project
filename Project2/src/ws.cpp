#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<sstream>
#include "armadillo"
#include "eig_solver.hpp"
#include "eig_output.hpp"
using namespace std;
using namespace arma;

const double pi=3.14159265359;
const double epsilon=1e-5; 
const int num=4; //only first num eigenvalues and corresponding eigenvectors are output
const double hbarc=197.3269788; 

int main(int argc, char* argv[])
{
    // command line input
    string filename; 
    int n; double rmax; 
    if (argc != 4)
    {
        cout << "Error in the input arguments! Please input three arguments: the number of points n, maximum radius rmax, the output file name."<<endl;
        return 0; 
    }
    else
    {
        n=atoi(argv[1]);
        if (n<=3)
        {
            cout << "Error! Too few points!"<<endl;
            return 0;
        }
        rmax=atof(argv[2]); 
        filename=argv[3];
    }
    
    int ia,iz; 
    cout <<"Input ia and iz:"<<endl; 
    cin >> ia >>iz; 
    int l,j2,tz; 
    cout <<"Input l and j*2:" <<endl; 
    cin >>l>>j2; 
    cout <<"Proton or neutron? (0 for neutron and 1 for proton)"<<endl; 
    cin >> tz; 
    if (tz>0) tz=1; else tz=0; 
    
    //initialize matrix
    double d,a,h,rho,pot; 
    h=rmax/n; 
    d=2/h/h;a=-1/h/h;  
    
    mat Amat(n-1,n-1);
    Amat.zeros();
    double v0,v1,r0,a0,vs,rc,rr0,rrc,esqr; 
    double mass; 
    v0=-51.0; v1=-33.0; r0=1.27; a0=0.67; vs=22.0; rc=1.2; 
    rr0=r0*pow(ia,1.0/3); rrc=rc*pow(ia,1.0/3); 
    esqr=hbarc/137; 
    if (tz>0) mass=1.00728*931.50*(ia-1)/ia; else mass=1.00866*931.50*(ia-1)/ia; 
    int fl; 
    fl=-(l+1);
    if (j2>2*l) fl=l;
    if (j2==2*l) fl=0.0;
    for (int i=0;i<n-1;i++)
    {
        rho=(i+1)*h; 
        Amat(i,i)=d+l*(l+1)/rho/rho; 
        pot=v0/(1.0+exp((rho-rr0)/a0))-vs*fl*exp((rho-rr0)/a0)/(a0*rho*pow(1+exp((rho-rr0)/a0),2));
        if (rho>rrc) 
            pot=pot+tz*iz*esqr/rho; 
        else
            pot=pot+tz*iz*esqr/rrc*(3.0/2.0-rho*rho/2.0/rrc/rrc); 
        Amat(i,i)=Amat(i,i)+pot/hbarc/hbarc*2*mass; 
        if (i+1<n-1) Amat(i,i+1)=a;
        if (i-1>=0) Amat(i,i-1)=a;
    }
    
    //armadillo eigenvalue solver
    vec eigenval(n-1); 
    mat eigenvec(n-1,n-1); 
    double start,finish; 
    start=clock();
    gen_mat_eig(n,Amat,eigenval,eigenvec);
    finish=clock();
    eigenval=eigenval*hbarc*hbarc/2/mass; 
    //eigenvec=eigenvec/sqrt(h); 
    output_all(filename+"_ws.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax,eigenval,eigenvec,num); 
    
    double rms; 
    for (int index=0;index<num;index++)
    {
        rms=0.0; 
        for (int i=0;i<n-1;i++)
        {
            rho=(i+1)*h; 
            rms=rms+eigenvec(i,index)*eigenvec(i,index)*rho*rho; 
        }
        rms=sqrt(rms); 
        cout <<"Energy: "<<eigenval(index)<<" rms: "<<rms<<endl; 
    }
    
    return 0;
}
