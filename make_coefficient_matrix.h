#ifndef MAKE_MATRIX
#define MAKE_MATRIX
#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
using T = Eigen::Triplet <double> ;
using sparse_matrix = Eigen::SparseMatrix <double> ;
using Solver = Eigen::SparseLU < sparse_matrix, Eigen::COLAMDOrdering<int> >;
void Make_Matrix(Solver **solver
                 ,sparse_matrix ***RHS_Coef_Matrix);//係数行列生成

double inner_product(int diff_flag,int half, double Delta,
                     int val, int test_order, int basis_order);
/*基底関数同士または基底関数とその導関数の内積を計算する
  関数を引数diff_flagに従って呼び出す*/

double element(
               int LHS_E_or_H,int  LHS_xyz, //未知数がHx,Hy,Ezのいずれであるか
               int  RHS_E_or_H,int  RHS_xyz, //Hx,Hy,Ezのいずれが係数にかかるか
               int *half, //境界条件によって基底関数が半分になるか否か
               int test_x_order,int  test_y_order, //試験関数の対応する導関数の次数
               int basis_x_order,int  basis_y_order, //基底関数の対応する導関数の次数
               int ix ,int iy, //試験関数の中心ノードの座標
               int x_val,int y_val //基底関数の中心ノードの座標 - 試験関数の中心ノードの座標
               );
/*
係数行列の各要素の値を返す。
*/
#endif



