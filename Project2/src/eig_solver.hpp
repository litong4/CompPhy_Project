#ifndef   EIG_SOLVER
#define   EIG_SOLVER

using namespace std; 
using namespace arma; 

bool jacobi_check(const mat&,const int,const double,int &,int &);
void gen_mat_eig(const int,const mat&,vec&,mat&);
bool gen_mat_jacobi(const int,mat,vec&,mat&,const double,const int);

#endif