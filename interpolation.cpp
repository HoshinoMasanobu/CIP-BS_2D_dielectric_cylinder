#include "set_variables.h"
#include "inner_product.h"
#include "make_coefficient_matrix.h"
#include "set_array.h"

double basis_func_cent0(double Delta, int node, int order,double x){
  double sign = 1.0;
  if(x < 0)sign = -1.0;
  double result = 0.0;
  for(int coefficient_order = K+1; coefficient_order <= 2*K+1; coefficient_order++){
    result += cf_a(order, coefficient_order )*std::pow(sign*x,coefficient_order);
  }
  result += std::pow(sign*x , order);
  result *= std::pow(Delta,order)/factorial(order);
  return  std::pow(sign,order)*result;
}

double EM_interpolated(Eigen::VectorXd **EM_density,int xyz,
                       double x, double y){
  int
    offset_ix,
    offset_iy,
    ix = 1, iy = 1;
  double
    offset_lx,
    offset_ly,
    DeltaX,
    DeltaY;
  // std::cout<<x <<" "<<y <<std::endl;
  if(x < RX[0] ){
    
    DeltaX = DX[0];
    offset_ix = 0;
    offset_lx = 0.0;
  }else if(x < RX[0] + RX[1]){
    DeltaX = DX[1];
    offset_ix = NX[0];
    offset_lx = RX[0];
  }else{
   
    DeltaX = DX[2];
    offset_ix = NX[1];
    offset_lx = RX[0] + RX[1];
  }
   
  if(y <  RY[0] ){
    
    DeltaY = DY[0];
    offset_iy = 0;
    offset_ly = 0.0;
  }else if(y < RY[0] + RY[1]){
    
    DeltaY = DY[1];
    offset_iy = NY[0];
    offset_ly = RY[0];
  }else{
    
    DeltaY = DY[2];
    offset_iy = NY[1];
    offset_ly = RY[0] + RY[1];
  }
   
  while(x > ix*DeltaX + offset_lx){
    ix++;
  }
  while(y > iy*DeltaY + offset_ly){
    iy++;
  }
  double sum = 0.0;
  
   for(int iterator = 0; iterator <= 3; iterator ++){
  
    int
      xx =  ix - iterator%2,
      yy =  iy - (iterator/2)%2;
    // std::cout<<ix<<" "<<iy<<" "<<xx<<" "<<yy<<std::endl;
    for(int bx = 0; bx <= K; bx++){
      for(int by = 0; by <= K; by ++){
        int node = Node( xyz, xx + offset_ix, yy + offset_iy, bx, by);
        //std::cout<<node<<std::endl;
        double EM = 0.0;
        
        if(node != -1)
          EM = (*EM_density)[xyz](node);
        
        sum +=
          EM
          *basis_func_cent0(DeltaX, xx, bx , (x -  (xx*DeltaX + offset_lx))/DeltaX)
          *basis_func_cent0(DeltaY, yy, by , (y -  (yy*DeltaY + offset_ly))/DeltaY);
        
      }
    }
  }
  return sum;
}



















