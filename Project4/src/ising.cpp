#include <iostream>
#include <random>
#include <cmath> 
#include <ctime> 
#include <cstdlib> 
#include <fstream>

using namespace std; 

int main(int argc, char **argv)
{
    int n,mc; 
    double temperature; 
    string filename; 
    //input
    if (argc != 5)
    {
        cout << "Error in the input arguments! Please input three arguments: size of the system N, temperature T, number of Monte Carlo cycles MC and output filename. "<<endl;
        return 0; 
    }
    else
    {
        n=atof(argv[1]);
        temperature=atof(argv[2]); 
        mc=atof(argv[3]); 
        filename=argv[4];
    }
    ofstream outfile; 
    outfile.open(filename+".txt"); 
    if (!outfile) 
    {
        cerr << "Cannot open output file!"; 
        return 1; 
    }
    
    //initialization
    int n_a=n+2; 
    enum spin:int {up=1, down=-1}; 
    int magnetic, energy; 
    magnetic=energy=0; 
    spin *a=new spin[n_a*n_a]; 
    srand(time(0)); 
    for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
        {
            a[i*n_a+j]=(double(rand())/RAND_MAX>0.5)?up:down; 
            magnetic+=a[i*n_a+j]; 
        }
    for (int i=1; i<=n; i++) //periodic boundary condition 
    {
        a[i]=a[n*n_a+i]; 
        a[(n+1)*n_a+i]=a[n_a+i]; 
        a[i*n_a]=a[i*n_a+n]; 
        a[i*n_a+(n+1)]=a[i*n_a+1]; 
    }
    for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
            energy-=a[i*n_a+j]*(a[i*n_a+(j+1)]+a[(i+1)*n_a+j]); 
    
    //Metropolis algorithm 
    double mag_avg,energy_avg,mag_tot,energy_tot; 
    mag_tot=magnetic; energy_tot=energy; 
    int flip,delta_e,delta_m;
    double exponential[5]; 
    for (int k=0; k<5; k++)
        exponential[k]=exp(-4.0*(k-2)/temperature); 
    for (int k=1; k<=mc; k++)
    {
        for (int i=1; i<=n; i++)
            for (int j=1; j<=n; j++)
            {
                flip=-a[i*n_a+j];  
                delta_m=2*flip; 
                delta_e=-(a[(i-1)*n_a+j]+a[(i+1)*n_a+j]+a[i*n_a+(j-1)]+a[i*n_a+(j+1)])*flip; 
                //if (abs(exponential[delta_e/2+2]-exp(-delta_e*2/temperature))>1e-6) //for debug
                //cout <<delta_e<<' '<<exponential[delta_e/2+2]<<' '<<exp(-delta_e*2/temperature)<<endl; 
                if ((delta_e<0)||(double(rand())/RAND_MAX<exponential[delta_e/2+2]))
                {
                    a[i*n_a+j]=spin(flip); 
                    if (i==1) a[(n+1)*n_a+j]=spin(flip); 
                    if (i==n) a[j]=spin(flip); 
                    if (j==1) a[i*n_a+(n+1)]=spin(flip); 
                    if (j==n) a[i*n_a]=spin(flip); 
                    magnetic+=delta_m; energy+=delta_e*2; 
                }
                mag_tot+=magnetic; energy_tot+=energy; 
                outfile <<magnetic<<' '<<energy<<endl; 
            }
    }
    mag_avg=mag_tot/mc/n/n; 
    energy_avg=energy_tot/mc/n/n; 
    
    //test output
    for (int i=0; i<n_a; i++)
    {
        for (int j=0; j<n_a; j++)
            cout <<a[i*n_a+j]<<' '; 
        cout <<endl; 
    }
    cout <<"magnetic: "<<magnetic<<' '<<mag_avg<<endl; 
    cout <<"energy: "<<energy<<' '<<energy_avg<<endl; 
    
    delete []a; 
    return 0; 
}
