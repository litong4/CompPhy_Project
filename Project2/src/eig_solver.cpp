#include <cmath>
#include "armadillo"
#include "eig_solver.hpp"

using namespace std; 
using namespace arma; 

//Armadillo eigenvalue decomposition
void gen_mat_eig(const int n,const mat& A,vec& eigval,mat& eigvec)
{
    eig_sym(eigval,eigvec,A); 
}

//check whether Jacobi algorithm should stop. If not, return the largest non-diagonal element's index (k,l) with k>l. 
bool jacobi_check(const mat& A,const int n,const double epsilon,int& k,int& l) //return k>l
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
bool gen_mat_jacobi(const int n,mat A,vec& eigval,mat& eigvec,const double epsilon,const int maxiteration)
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
