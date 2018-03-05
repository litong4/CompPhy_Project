#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "armadillo"
#include "eig_solver.hpp"
#include "eig_output.hpp"
using namespace std;
using namespace arma;

const double pi=3.14159265359;
const double epsilon=1e-5; 
int maxiteration=1000000; 

int main(int argc, char* argv[])
{
    // command line input
    string filename; 
    int n;
    if (argc != 3)
    {
        cout << "Error in the input arguments! Please input two arguments: the number of points n and the output file name."<<endl;
        return 0; 
    }
    else
    {
        n=atoi(argv[1]);
        if (n<=2)
        {
            cout << "Error! Too few points!"<<endl;
            return 0;
        }
        filename=argv[2];
    }
    
    //initialize matrix elements
    double d,a; 
    d=2*n*n;a=-1*n*n;  
    
    //initialization
    mat Amat(n-1,n-1);
    Amat.zeros();
    for (int i=0;i<n-1;i++)
    {
        Amat(i,i)=d;
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
    output_all(filename+"_arma.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,1.0,eigenval,eigenvec,-1); 
    
    //Jacobi's method
    start=clock();
    if (gen_mat_jacobi(n,Amat,eigenval,eigenvec,epsilon,maxiteration))
    {
        finish=clock();
        output_all(filename+"_jacobi.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,1.0,eigenval,eigenvec,-1); 
    }
    else
    {
        finish=clock(); 
        output_fail(filename+"_jacobi.txt","Jacobi method exceeds maximum iteration number!",(double)(finish-start)/CLOCKS_PER_SEC,n,1.0); 
    }
    
    return 0;
}
