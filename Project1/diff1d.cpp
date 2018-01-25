#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>
using namespace std;

double fun_f(double);
double fun_ans(double);
void gen_tri_solve(int,const double*,const double*,const double*, double*, const double*);
void output_all(const string&,double,int,const double*,const double*,const double*);
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
        filename=argv[2];
    }
    
    //construct vectors and initialization
    double h = 1.0/n;    
    double *x, *y; // for correct answers
    double *u, *b; // for Au=b
    double *d, *e, *f; // for (n-1)*(n-1) tri-diagonal matrix A, d is diagonal, e is lower, f is upper
    x=new double[n+1]; y=new double[n+1];
    u=new double[n+1]; b=new double[n+1];
    d=new double[n+1]; e=new double[n+1]; f=new double[n+1]; 
        //d[0], d[n], e[0], e[n], e[n-1], f[0], f[n], f[n-1] are never used
    u[0]=0.0; u[n]=0.0; 
    for (int i=0; i<=n; i++)
    {
        x[i]=h*i; 
        b[i]=h*h*fun_f(x[i]); 
        y[i]=fun_ans(x[i]);
        d[i]=2; e[i]=-1; f[i]=-1; 
    }
    
    clock_t start, finish;
    //use a general solver for tri-diagonal matrix 
    start=clock(); 
    gen_tri_solve(n,d,e,f,u,b);
    finish=clock(); 
    output_all(filename,(double)(finish-start)/CLOCKS_PER_SEC,n,x,u,y);
    
    //delete
    delete[] x; delete[] y; delete[] u; delete[] b;
    delete[] d; delete[] e; delete[] f; 
    
    return 0;
}

inline double fun_f(double x)
{
    return (100*exp(-10*x)); 
}

inline double fun_ans(double x)
{
    return (1-(1-exp(-10))*x-exp(-10*x)); 
}

void output_all(const string& filename,double time,int n,const double* x,const double* u,const double* y)
{
    ofstream outfile;
    double rel_error, avg_error=0.0;
    
    outfile.open(filename);
    outfile <<"n= "<<n<<endl;
    outfile <<"Use time "<<time<<" second."<<endl;
    outfile <<"x, u, y, relative error"<<endl;
    outfile <<x[0]<<' '<<u[0]<<' '<<y[0]<<' '<<0.0<<endl;
    for (int i=1;i<n;i++)
    {
        rel_error=abs((u[i]-y[i])/y[i]); 
        outfile <<x[i]<<' '<<u[i]<<' '<<y[i]<<' '<<rel_error<<endl;
        avg_error=avg_error+rel_error;
    }
    outfile <<x[n]<<' '<<u[n]<<' '<<y[n]<<' '<<0.0<<endl;
    avg_error=avg_error/(n-1);
    outfile <<"Average relative error is "<<avg_error<<endl; 
    outfile <<endl; 
    outfile.close(); 
}

void gen_tri_solve(int n,const double* d,const double* e,const double* f, double* u, const double* b)
{
    double *d_t, *b_t; 
    double temp;
    d_t=new double[n+1]; b_t=new double[n+1];
    
    //forward substitution
    d_t[1]=d[1]; b_t[1]=b[1]; 
    for (int i=2; i<n; i++)
    {
        temp=e[i-1]/d_t[i-1];
        d_t[i]=d[i]-f[i-1]*temp;
        b_t[i]=b[i]-b_t[i-1]*temp;
    }
    
    //backward substitution
    u[n-1]=b_t[n-1]/d_t[n-1];
    for (int i=n-2; i>0; i--)
    {
        u[i]=(b_t[i]-e[i]*u[i+1])/d_t[i];
    }
    
    delete[] d_t; delete[] b_t;
}