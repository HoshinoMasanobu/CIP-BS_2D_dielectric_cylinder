#ifndef ARRAY
#define ARRAY
#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
extern int
*****node;//[E_or_H (0,1)][x y or z (0,1,2)][x][y][z][x_order][y_order][z_order]
extern int
*row_size;//[E_or_H][x_y_or_z] ex)Ez = [0][2],Hy = [1][1]
int relative_index(int i,int *ii/*minus i+1*/,int *iii/*plus i-1*/);
void set_node_number();
int Node(int x_y_or_z, int x, int y, int x_order,int y_order);
int Boundary_condition(int E_or_H,
                       int x_y_z,
                       int x,int y,
                       int x_order,int y_order);
#endif 
