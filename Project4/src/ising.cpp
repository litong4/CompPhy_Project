#include <iostream>
#include <random>
#include <cmath> 
#include <ctime> 
#include <cstdlib> 
#include <fstream>

using namespace std; 

int main(int argc, char **argv)
{
    int n,mc,ini,out_detail; 
    double temperature; 
    string filename; 
    //input
    if (argc != 7)
    {
        cout << "Error in the input arguments! \nPlease input five arguments: size of the system N, temperature T, number of Monte Carlo cycles MC, whether random initialization (0 or 1), whether output Monte Carlo process (0 or 1) and output filename. "<<endl;
        return 0; 
    }
    else
    {
        n=atof(argv[1]);
        temperature=atof(argv[2]); 
        mc=atof(argv[3]); 
        ini=atof(argv[4]); 
        out_detail=atof(argv[5]); 
        filename=argv[6];
    }
    ofstream outfile,summary; 
    if (out_detail) 
    {
        outfile.open(filename+"_mc.txt"); 
        if (!outfile)
        {
            cerr << "Cannot open output file!"; 
            return 1; 
        }
        outfile <<"Magnetization M and energy E per spin in the Monte Carlo process"<<endl; 
    }
    summary.open(filename+"_sum.txt"); 
    if (!summary) 
    {
        cerr << "Cannot open output file!"; 
        return 1; 
    }
    summary <<"Summary of Monte-Carlo mean values of Ising model"<<endl; 
    
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
            if (ini) 
                a[i*n_a+j]=(double(rand())/RAND_MAX>0.5)?up:down; //random initialization
            else
                a[i*n_a+j]=up; //all point up
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
    double mag_avg,energy_avg,mag_sqr_avg,energy_sqr_avg; 
    double mag_tot,energy_tot,mag_sqr_tot,energy_sqr_tot; 
    mag_tot=0;energy_tot=0;mag_sqr_tot=0;energy_sqr_tot=0; 
    int flip,delta_e,delta_m;
    double exponential[5]; 
    for (int k=0; k<5; k++)
        exponential[k]=exp(-4.0*(k-2)/temperature); 
    
    random_device rd; 
    mt19937 gen(rd()); 
    double gen_max=gen.max(); 
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
                if ((delta_e<0)||(gen()/gen_max<exponential[delta_e/2+2]))
                //if ((delta_e<0)||(gen()/gen_max<exp(-delta_e*2/temperature)))
                {
                    a[i*n_a+j]=spin(flip); 
                    if (i==1) a[(n+1)*n_a+j]=spin(flip); 
                    if (i==n) a[j]=spin(flip); 
                    if (j==1) a[i*n_a+(n+1)]=spin(flip); 
                    if (j==n) a[i*n_a]=spin(flip); 
                    magnetic+=delta_m; energy+=delta_e*2; 
                }
                //if (mc>mc/10) 
                //{
                mag_tot+=abs(magnetic); energy_tot+=energy; 
                mag_sqr_tot+=magnetic*magnetic; energy_sqr_tot+=energy*energy; 
                //}
                if (outfile) outfile <<double(magnetic)/n/n<<' '<<double(energy)/n/n<<endl; 
            }
    }
    //mc=mc-mc/10; 
    mag_avg=mag_tot/mc/n/n; 
    energy_avg=energy_tot/mc/n/n; 
    mag_sqr_avg=mag_sqr_tot/mc/n/n; 
    energy_sqr_avg=energy_sqr_tot/mc/n/n; 
    
    double chi,cv; 
    cv=(energy_sqr_avg-energy_avg*energy_avg)/temperature/temperature; 
    chi=(mag_sqr_avg-mag_avg*mag_avg)*temperature; 
    
    summary <<"Magnetization: "<<mag_avg<<endl; 
    summary <<"Magnetization^2: "<<mag_sqr_avg<<endl; 
    summary <<"Energy: "<<energy_avg<<endl; 
    summary <<"Energy^2: "<<energy_sqr_avg<<endl; 
    summary <<"Specific heat: "<<cv<<endl; 
    summary <<"Susceptibility: "<<chi<<endl; 
    
    if (n==2) 
    {
        double cch,ssh,eep; 
        cch=cosh(8.0/temperature); ssh=sinh(8.0/temperature); eep=exp(8.0/temperature); 
        double energy_th,mag_th,cv_th,energy_sqr_th,mag_sqr_th,chi_th; 
        mag_th=(2.0*eep+1.0)/(cch+3.0); 
        mag_sqr_th=8.0*(eep+1.0)/(cch+3.0); 
        energy_th=-8.0*ssh/(cch+3); 
        energy_sqr_th=64.0*cch/(3.0+cch); 
        cv_th=64*(1+3.0*cch)/(3.0+cch)/(3.0+cch)/temperature/temperature; 
        chi_th=(mag_sqr_th-mag_th*mag_th)*temperature; 
        summary <<"Magnetization (analytical): "<<mag_th<<endl; 
        summary <<"Magnetization^2 (analytical): "<<mag_sqr_th<<endl; 
        summary <<"Energy (analytical): "<<energy_th<<endl; 
        summary <<"Energy^2 (analytical): "<<energy_sqr_th<<endl; 
        summary <<"Specific heat (analytical): "<<cv_th<<endl; 
        summary <<"Susceptibility (analytical): "<<chi_th<<endl; 
    }
    
    delete []a; 
    if (outfile) outfile.close(); 
    summary.close(); 
    
    return 0; 
}
