#include <iostream>
#include "armadillo"
using namespace std;

int main()
{
arma::arma_version ver;
std::cout << "ARMA version: "<< ver.as_string() << std::endl;
}
