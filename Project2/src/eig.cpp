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
const double epsilon=1e-5; 
int maxiteration=1000000; 

void output_all(const string&,const double,const int,const double,const double,vec&,mat&);
void output_fail(const string&,const string&,const double,const int,const double,const double);
void gen_mat_eig(const int,const mat&,vec&,mat&);
bool gen_mat_jacobi(const int,mat,vec&,mat&);

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
    if (gen_mat_jacobi(n,Amat,eigenval,eigenvec))
    {
        finish=clock();
        output_all(filename+"_jacobi.txt",(double)(finish-start)/CLOCKS_PER_SEC,n,d,a,eigenval,eigenvec); 
    }
    else
    {
        finish=clock(); 
        output_fail(filename+"_jacobi.txt","Jacobi method exceeds maximum iteration number!",(double)(finish-start)/CLOCKS_PER_SEC,n,d,a); 
    }
    
    return 0;
}

//Armadillo eigenvalue decomposition
void gen_mat_eig(const int n,const mat& A,vec& eigval,mat& eigvec)
{
    eig_sym(eigval,eigvec,A); 
}

//check whether Jacobi algorithm should stop. If not, return the largest non-diagonal element's index (k,l) with k>l. 
bool jacobi_check(const mat& A,const int n,const double epsilon,int& k,int& l)
{
    double maximum=-1.0; 
    for (int i=0;i<n-1;i++)
        for (int j=0;j<i;j++)
        {
            if (abs(A(i,j))>maximum) 
            {
                maximum=abs(A(i,j)); 
                k=i;l=j; 
            }
        }
    
    if (maximum>epsilon) 
        return 1;
    else
        return 0; 
}

//Jacobi algorithm for eigenvalue decomposition 
bool gen_mat_jacobi(const int n,mat A,vec& eigval,mat& eigvec)
{
    int k,l; 
    double tau,t,c,s,temp_ik,temp_il,a_kk,a_ll,a_kl;
    int iter=0; 
    bool flag=true; 
    
    eigvec.eye(); 
    while (jacobi_check(A,n,epsilon,k,l))   //k>l
    {
        iter++; 
        if (iter>maxiteration)
        {
            flag=false;
            break;
        }
        
        tau=(A(l,l)-A(k,k))/2.0/A(k,l); 
        if (tau>=0)
            t=1.0/(tau+sqrt(1.0+tau*tau));
        else
            t=-1.0/(-tau+sqrt(1.0+tau*tau)); 
        c=1.0/sqrt(1+t*t); 
        s=t*c; 

        a_kk=A(k,k); a_ll=A(l,l); a_kl=A(k,l); 
        A(k,k)=a_kk*c*c-2*a_kl*c*s+a_ll*s*s;
        A(l,l)=a_ll*c*c+2*a_kl*c*s+a_kk*s*s; 
        A(l,k)=0.0; A(k,l)=0.0; 
        for (int i=0;i<n-1;i++)
        {
            if ((i!=k)&&(i!=l))
            {
                temp_ik=A(i,k); temp_il=A(i,l); 
                A(i,k)=temp_ik*c-temp_il*s; 
                A(i,l)=temp_il*c+temp_ik*s; 
                A(k,i)=A(i,k); A(l,i)=A(i,l); 
            }
            temp_il=eigvec(i,l); temp_ik=eigvec(i,k); 
            eigvec(i,k)=-temp_il*s+temp_ik*c; 
            eigvec(i,l)=temp_il*c+temp_ik*s; 
        }
    }
    
    for (int i=0;i<n-1;i++)
        eigval(i)=A(i,i); 
    
    return flag;
}

void output_all(const string& filename,const double time,const int n,const double d,const double a,vec& eigval,mat& eigvec)
{
    ofstream outfile;
    double rel_error, max_error=-1.0;
    double eigval_correct; 
    uvec index_order(n-1); 
    
    outfile.open(filename);
    outfile <<"n = "<<n<<endl;
    outfile <<"d = "<<d<<" and a = "<<a<<endl;
    outfile <<"Use time "<<time<<" seconds."<<endl; 
    
    index_order=sort_index(eigval); 
    eigval=sort(eigval);  
    eigvec=eigvec.cols(index_order); 
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

void output_fail(const string& filename,const string& message,const double time,const int n,const double d,const double a)
{
    ofstream outfile;
    outfile.open(filename);    
    outfile <<message<<endl; 
    outfile <<"n = "<<n<<endl;
    outfile <<"d = "<<d<<" and a = "<<a<<endl;
    outfile <<"Use time "<<time<<" seconds."<<endl; 
}
