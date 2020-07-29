#include <Eigen/Sparse>
#include "set_variables.h"
#include "initialize.h"
#include "set_array.h"
#include "gaussian_pulse.h"
#include "ununiform.h"
void initialize(Eigen::VectorXd **EM_density){
  ////////////////////初期値入力///////////////////
  
  for(int x = 0; x <= NX[2]; x++){
    for(int y = 0; y <= NY[2]; y++){
      
      for(int x_order = 0; x_order <= K;  x_order ++){//m
        for(int y_order = 0; y_order <= K;  y_order ++){
          
          
          for(int xyz = 0; xyz <= 2; xyz++){
            int node_num =
              Node( xyz,  x,  y, x_order, y_order);
            if(node_num != -1){
              (*EM_density)[xyz](node_num) = 0.0;
              if(xyz == 2)
                (*EM_density)[xyz](node_num) = 0.0;
                }
          }
        }
      }
    }
  }
} 

















