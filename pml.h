#ifndef PML
#define PML
#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
extern
int
  **pml_start,
  **pml_len;

double sigma(int xy, int i, int diff);
void set_node_num_pml();//pmlの境界のノード番号格納
void update_H_pml(Eigen::VectorXd **EM_density,Eigen::VectorXd *EM_past,
                  Eigen::VectorXd **pmlH);//補助場pmlH計算
void update_E_pml(Eigen::VectorXd **EM_density,Eigen::VectorXd *EM_past,
                  Eigen::VectorXd **pmlEz);//補助場pmlEz計算
#endif



