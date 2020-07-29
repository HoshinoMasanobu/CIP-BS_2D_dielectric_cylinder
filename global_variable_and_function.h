#ifndef HEADER
#define HEADER
#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
using T = Eigen::Triplet <double> ;
using sparse_matrix = Eigen::SparseMatrix <double> ;
using Solver = Eigen::SparseLU < sparse_matrix, Eigen::COLAMDOrdering<int> >;

int node_number(int,int,int,int,int);//time, position, time_order, x_order, E or H
double cf_a(int,int);
void set_variables(void);
void set_node_number(void);
int relative_index(int i,int *ii,int *iii);
int Node(int E_or_H,int x_y_or_z, int x, int y,int z,int x_order,int y_order,int z_order);
void Make_Matrix(Solver ***solver
                 ,sparse_matrix ****RHS_Coef_Matrix);
double factorial(int);
double combination(int,int);
void basis_set(void);
double inner_product_S(double,int,int,int);
double inner_product_L(double,int,int,int);
double Gaussian_Pulse(double eta,double xi,double sigma,int order);
void initialize(Eigen::VectorXd ***EM_density);
void free_a(void);

extern int
**row_size;//[E_or_H][x_y_or_z] ex)Ez = [0][2],Hy = [1][1]
extern int a_matrix_size;
constexpr int K { 1 };
constexpr double 
c { 3.0e8 },//光速
  mu { 4.0 * M_PI *1.0e-7 },//μ、透磁率
  eps { 1.0 / (mu * c * c) };//ε、誘電率
extern int
NT,
  NX,
  NY,
  NZ,
  Sorce_point_x,
  Sorce_point_y,
  Sorce_point_z;
extern double
  C,//クーラン数 
  RX,//領域の長さ x
  RY,//y
  RZ,//z
  DX,
  DY,
  DZ,
  DT;
extern Eigen::VectorXd* a;
extern int
********node;//[E_or_H (0,1)][x y or z (0,1,2)][x][y][z][x_order][y_order][z_order]
#endif

