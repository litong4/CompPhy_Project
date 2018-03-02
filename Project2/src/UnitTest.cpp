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
    const double d=3, a=-2, pi=3.14159265359; 
    double fnorm_squ=0.0;  
    double temp; 
    mat Amat(n-1,n-1); 
    
    //intialization
    Amat.zeros(); 
    for (int i=0;i<n-1;i++)
    {
        //fnorm_squ is the square of F norm before diagonalization
        Amat(i,i)=d; 
        fnorm_squ=fnorm_squ+d*d; 
        if (i+1<n-1) {Amat(i,i+1)=a; fnorm_squ=fnorm_squ+a*a; }
        if (i-1>=0) {Amat(i,i-1)=a; fnorm_squ=fnorm_squ+a*a; }
    }
     
    vec eigval(n-1); mat eigvec(n-1,n-1);  
    gen_mat_eig(n,Amat,eigval,eigvec); //test armadillo eigenvalue solver
    
    temp=0.0; 
    for (int i=0;i<n-1;i++)
    {
        REQUIRE(eigval(i)==Approx(d+2*a*cos((i+1)*pi/n))); //test eigenvalues
        temp=temp+eigval(i)*eigval(i); //calculate square of F norm after diagonalization
    }
    REQUIRE(temp==Approx(fnorm_squ)); //test conservation of F norm
    
    //test eigenvectors
    mat testMat; 
    //first test orthonormalization of eigenvectors
    testMat=eigvec.t()*eigvec; 
    for (int i=0;i<n-1;i++)
        for (int j=0;j<n-1;j++)
            if (i==j) 
                REQUIRE(testMat(i,i)==Approx(1)); 
            else
                REQUIRE(abs(testMat(i,j))<1e-8);
    //then test the exact values of eigenvectors
    vec testvec(n-1); 
    testvec(0)=-0.5; testvec(1)=0.7071; testvec(2)=-0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(2)))==Approx(1)); 
    testvec(0)=-0.7071; testvec(1)=0.0; testvec(2)=0.7071;
    REQUIRE(abs(dot(testvec,eigvec.col(1)))==Approx(1)); 
    testvec(0)=0.5; testvec(1)=0.7071; testvec(2)=0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(0)))==Approx(1)); 
    
    double epsilon=1e-5; 
    int maxiteration=1000;
    REQUIRE(gen_mat_jacobi(n,Amat,eigval,eigvec,epsilon,maxiteration)==1); //test Jacobi eigenvalue solver
    uvec index_order=sort_index(eigval); 
    eigval=sort(eigval);  
    eigvec=eigvec.cols(index_order);
    
    temp=0.0; 
    for (int i=0;i<n-1;i++)
    {
        REQUIRE(eigval(i)==Approx(d+2*a*cos((i+1)*pi/n))); //test eigenvalues
        temp=temp+eigval(i)*eigval(i); //calculate square of F norm after diagonalization
    }
    REQUIRE(temp==Approx(fnorm_squ)); //test conservation of F norm
    
    //first test orthonormalization of eigenvectors
    testMat=eigvec.t()*eigvec; 
    for (int i=0;i<n-1;i++)
        for (int j=0;j<n-1;j++)
            if (i==j) 
                REQUIRE(testMat(i,i)==Approx(1)); 
            else
                REQUIRE(abs(testMat(i,j))<1e-8);
    //then test the exact values of eigenvectors
    testvec(0)=-0.5; testvec(1)=0.7071; testvec(2)=-0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(2)))==Approx(1)); 
    testvec(0)=-0.7071; testvec(1)=0.0; testvec(2)=0.7071;
    REQUIRE(abs(dot(testvec,eigvec.col(1)))==Approx(1)); 
    testvec(0)=0.5; testvec(1)=0.7071; testvec(2)=0.5;
    REQUIRE(abs(dot(testvec,eigvec.col(0)))==Approx(1)); 
}