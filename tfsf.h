#ifndef TFSF
#define TFSF

#include "make_coefficient_matrix.h"
extern double omega,phi_i;
void set_TFSF(sparse_matrix ***Coef_Mat);//入射波に乗ずる係数行列生成
void incident_update( int xyz, int n, Eigen::VectorXd *Incident);//入射波更新
double plane_wave(int xyz, int n, int ix ,int iy, int x_order,int y_order);//Gaussian pulse 使わない.



#endif
