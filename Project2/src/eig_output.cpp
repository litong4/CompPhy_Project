#include <fstream>
#include "armadillo"
#include "eig_output.hpp"

using namespace std; 
using namespace arma; 

//output the result of eigenvalue decomposition when the decomposition is done successfully. 
//parameter "howmany" gives how many smallest eigenvalues and corresponding eigenvectors should be output. 
//howmany<=0 stands for output all the results. 
void output_all(const string& filename,const double time,const int n,const double d,const double a,vec& eigval,mat& eigvec,int howmany)
{
    ofstream outfile;
    double rel_error, max_error=-1.0;
    double eigval_correct; 
    uvec index_order(n-1); 
    
    outfile.open(filename);
    outfile <<"n = "<<n<<endl;
    outfile <<"Use time "<<time<<" seconds."<<endl; 
    
    index_order=sort_index(eigval); 
    eigval=sort(eigval);  
    eigvec=eigvec.cols(index_order); 
    if (howmany>n-1) howmany=-1; 
    if (howmany<=0)
    {
        eigval.print(outfile,"Eigenvalues: ");
        eigvec.print(outfile,"Eigenvectors: "); 
    }
    else
    {
        eigval.subvec(1,howmany).print(outfile,"Eigenvalues: ");
        eigvec.cols(1,howmany).print(outfile,"Eigenvectors: ");         
    }
    
    outfile.close();
}

//output failure message when eigenvalue decomposition fails. 
void output_fail(const string& filename,const string& message,const double time,const int n,const double d,const double a)
{
    ofstream outfile;
    outfile.open(filename);    
    outfile <<message<<endl; 
    outfile <<"n = "<<n<<endl;
    outfile <<"Use time "<<time<<" seconds."<<endl; 
}