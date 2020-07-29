#ifndef INTERPOLATED
#define INTERPOLATED
#include <iostream>
#include <Eigen/Sparse>

double EM_interpolated(Eigen::VectorXd **EM_density, int xyz,
                       double x, double y);
double basis_func_cent0(double Delta, int node, int order,double x);

#endif
