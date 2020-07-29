#ifndef INNER_PRODUCT
#define INNER_PRODUCT

#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
double factorial(int n);

//2K+1次の基底関数の各係数
double cf_a(int basis_function_order,int coefficient_order);

//基底関数と基底関数の内積
double inner_product_S(int half, //試験関数の中心ノードから見て完全導体が-側にあるなら-1,+側にあるなら1、無いなら0
                       int prev_next, //試験関数の中心ノードから見て0なら-側,1なら+側片側だけの積分
                       double Delta_prev,double Delta_next, //不均一な分割によるセルサイズ
                       int w, int m, int n //w:(基底関数の中心ノード)-(試験関数の中心ノード), m:試験関数の次数,n:基底関数の次数
                       );

//基底関数と基底関数の導関数の内積
double inner_product_L(int half,
                       int prev_next,
                       double Delta_prev,double Delta_next,
                       int w, int m, int n);

void basis_set();
#endif
