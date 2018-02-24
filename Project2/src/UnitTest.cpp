#define CATCH_CONFIG_MAIN 

#include <cmath>
#include "armadillo"
#include "catch.hpp"
#include "eig_solver.hpp"

using namespace std; 
using namespace arma; 

TEST_CASE("Jacobi check is performed to find the maximum off-diagonal element in a matrix","[Jacobi check]")
{
    const int n=5; 
    mat A(n-1,n-1); 
    A(0,0)=1.2; A(0,1)=2.1; A(0,2)=2.2; A(0,3)=-4.1; 
    A(1,0)=2.1; A(1,1)=-5.5; A(1,2)=0.3; A(1,3)=-2.2; 
    A(2,0)=2.2; A(2,1)=0.3; A(2,2)=4.5; A(2,3)=0.5; 
    A(3,0)=-4.1; A(3,1)=-2.2; A(3,2)=0.5; A(3,3)=7.0; 
    int k,l; 
    double epsilon=1e-5; 
    REQUIRE(jacobi_check(A,n,epsilon,k,l)==1); 
    REQUIRE(k>l); 
    REQUIRE(k==3); 
    REQUIRE(l==0); 
}

TEST_CASE("Eigenvalue and eigenvector test using a tridiagonal Toeplitz matrix","[Eigenvalue and eigenvector test]")
{
    const int n=4; 
    mat A(n-1,n-1); 
    A.zeros(); 
    for (int i=0;i<n-1;i++)
    {
        A(i,i)=3; 
        if (i+1<n-1) A(i,i+1)=-2; 
        if (i-1>=0) A(i,i-1)=-2; 
    }
     
    vec eigval(n-1); mat eigvec(n-1,n-1);  
    gen_mat_eig(n,A,eigval,eigvec); 
    REQUIRE(eigval(0)==Approx(0.17157).epsilon(1e-3)); 
    REQUIRE(eigval(1)==Approx(3.00000).epsilon(1e-3)); 
    REQUIRE(eigval(2)==Approx(5.82843).epsilon(1e-3)); 
    vec testvec(n-1); 
    testvec(0)=-0.5; testvec(1)=0.7071; testvec(2)=-0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(2)))==Approx(1)); 
    testvec(0)=-0.7071; testvec(1)=0.0; testvec(2)=0.7071;
    REQUIRE(abs(dot(testvec,eigvec.col(1)))==Approx(1)); 
    testvec(0)=0.5; testvec(1)=0.7071; testvec(2)=0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(0)))==Approx(1)); 
    
    double epsilon=1e-5; 
    int maxiteration=1000;
    REQUIRE(gen_mat_jacobi(n,A,eigval,eigvec,epsilon,maxiteration)==1); 
    uvec index_order=sort_index(eigval); 
    eigval=sort(eigval);  
    eigvec=eigvec.cols(index_order);
    REQUIRE(eigval(0)==Approx(0.17151).epsilon(1e-3)); 
    REQUIRE(eigval(1)==Approx(3.00000).epsilon(1e-3)); 
    REQUIRE(eigval(2)==Approx(5.82843).epsilon(1e-3));    
    testvec(0)=-0.5; testvec(1)=0.7071; testvec(2)=-0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(2)))==Approx(1)); 
    testvec(0)=-0.7071; testvec(1)=0.0; testvec(2)=0.7071;
    REQUIRE(abs(dot(testvec,eigvec.col(1)))==Approx(1)); 
    testvec(0)=0.5; testvec(1)=0.7071; testvec(2)=0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(0)))==Approx(1)); 
}