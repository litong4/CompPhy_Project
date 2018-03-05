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
int maxiteration=1000000; 
const int num=4; //only first num eigenvalues and corresponding eigenvectors are output

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
    
    //------------------------one-electron case------------------------//
    //initialize matrix
    double d,a,h; 
    h=rmax/n; 
    d=2/h/h;a=-1/h/h;  
    
    mat Amat(n-1,n-1);
    Amat.zeros();
    for (int i=0;i<n-1;i++)
    {
        Amat(i,i)=d+(i+1)*h*(i+1)*h;
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
    output_all(filename+"_one_arma.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax,eigenval,eigenvec,num); 
    
    //Jacobi's method
    start=clock();
    if (gen_mat_jacobi(n,Amat,eigenval,eigenvec,epsilon,maxiteration))
    {
        finish=clock();
        output_all(filename+"_one_jacobi.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax,eigenval,eigenvec,num); 
    }
    else
    {
        finish=clock(); 
        output_fail(filename+"_one_jacobi.txt","Jacobi method exceeds maximum iteration number!",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax); 
    }
    
    //------------------------two-electron case------------------------//
    
    double wr[4]; 
    wr[0]=0.01; wr[1]=0.5; wr[2]=1; wr[3]=5; 
    for (int index=0;index<4;index++)
    {
        //initialize matrix elements
        mat Amat(n-1,n-1);
        Amat.zeros();
        for (int i=0;i<n-1;i++)
        {
            Amat(i,i)=d+(i+1)*h*(i+1)*h*wr[index]*wr[index]+1.0/h/(i+1);
            if (i+1<n-1) Amat(i,i+1)=a;
            if (i-1>=0) Amat(i,i-1)=a;
        }
        
        stringstream outwr; 
        outwr << wr[index]; 
        
        //armadillo eigenvalue solver
        vec eigenval(n-1); 
        mat eigenvec(n-1,n-1); 
        double start,finish; 
        start=clock();
        gen_mat_eig(n,Amat,eigenval,eigenvec);
        finish=clock();
        output_all(filename+"_two_wr"+outwr.str()+"_arma.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax,eigenval,eigenvec,num); 
        
        //Jacobi's method
        start=clock();
        if (gen_mat_jacobi(n,Amat,eigenval,eigenvec,epsilon,maxiteration))
        {
            finish=clock();
            output_all(filename+"_two_wr"+outwr.str()+"_jacobi.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax,eigenval,eigenvec,num); 
        }
        else
        {
            finish=clock(); 
            output_fail(filename+"_two_wr"+outwr.str()+"_jacobi.txt","Jacobi method exceeds maximum iteration number!",(double)(finish-start)/CLOCKS_PER_SEC,n,rmax); 
        }
    }
    
    return 0;
}
