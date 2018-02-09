#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "armadillo"
using namespace std;
using namespace arma;

const double pi=3.14159265359;
const double epsilon=1e-6; 

void output_all(const string&,const double,const int,const double,const double,vec&,const mat&);
void gen_mat_eig(const int,const mat&,vec&,mat&);
void gen_mat_jacobi(const int,mat,vec&,mat&);

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
    
    //keyboard input
    double d,a;
    cout << "Please intput diagonal matrix element d and off-diagonal element a to make a tri-diagonal matrix:" << endl;
    cin >>d>>a;
    
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
    output_all(filename+"_arma.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,d,a,eigenval,eigenvec); 
    
    //Jacobi's method
    start=clock();
    gen_mat_jacobi(n,Amat,eigenval,eigenvec);
    finish=clock();
    eigenval=sort(eigenval); 
    output_all(filename+"_jacobi.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,d,a,eigenval,eigenvec); 
    
    return 0;
}

void gen_mat_eig(const int n,const mat& A,vec& eigval,mat& eigvec)
{
    eig_sym(eigval,eigvec,A); 
}

bool jacobi_check(const mat& A,const int n,const double epsilon,int& k,int& l)
{
    double norm=0.0; 
    double maximum=-1.0; 
    for (int i=0;i<n-1;i++)
        for (int j=0;j<i;j++)
        {
            norm=norm+A(i,j)*A(i,j)*2; 
            if (abs(A(i,j))>maximum) 
            {
                maximum=abs(A(i,j)); 
                k=i;l=j; 
            }
        }
    
    if (norm>epsilon) 
        return 1;
    else
        return 0; 
}

void gen_mat_jacobi(const int n,mat A,vec& eigval,mat& eigvec)
{
    int k,l; 
    double tau,t,c,s; 
    eigvec.eye(); 
    mat temp(n-1,n-1); 
    while (jacobi_check(A,n,epsilon,k,l))   //k>l
    {
        tau=(A(l,l)-A(k,k))/2.0/A(k,l); 
        t=-tau+sqrt(1.0+tau*tau); 
        c=1.0/sqrt(1+t*t); 
        s=t*c; 
        temp.eye(); 
        temp(k,k)=c; temp(l,l)=c; 
        temp(k,l)=-s; temp(l,k)=s; 
        eigvec=temp*eigvec; 
        A=temp.t()*A*temp; 
    }
    
    for (int i=0;i<n-1;i++)
        eigval(i)=A(i,i); 
}

void output_all(const string& filename,const double time,const int n,const double d,const double a,vec& eigval,const mat& eigvec)
{
    ofstream outfile;
    double rel_error, max_error=-1.0;
    double eigval_correct; 
    
    outfile.open(filename);
    outfile <<"n = "<<n<<endl;
    outfile <<"d = "<<d<<" and a = "<<a<<endl;
    outfile <<"Use time "<<time<<" seconds."<<endl; 
    
    eigval.print(outfile,"Eigenvalues: ");
    eigvec.print(outfile,"Eigenvectors: ");
    
    for (int i=0;i<n-1;i++)
    {
        eigval_correct=d+2*a*cos((i+1)*pi/n); 
        rel_error=abs((eigval(i)-eigval_correct)/eigval_correct); 
        if (rel_error>max_error) max_error=rel_error; 
        //outfile <<eigval_correct<<endl;
    }
    outfile<<"Maximum relative error is "<<max_error<<endl; 
    
    outfile.close();
}
