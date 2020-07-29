#include "make_coefficient_matrix.h"
#include "set_variables.h"
#include "set_array.h"
#include "inner_product.h"
#include "dielectrics.h"
#include "ununiform.h"
double inner_product(bool diff_flag,int half,
                     int prev_next,
                     double Delta_prev,double Delta_next,
                     int val, int test_order, int basis_order){
  
  if(diff_flag == 1)
    return inner_product_S(half,
                           prev_next,
                           Delta_prev,Delta_next,
                           val, test_order, basis_order);
  if(diff_flag == 0)
    return inner_product_L(half,
                           prev_next,
                           Delta_prev,Delta_next,
                           val, test_order, basis_order);
}

  std::string char_EH[2] = {"E","H"},char_xyz[3]={"x","y","z"};
bool diff_flag[3][3] =  {{1,0,0},{0,1,0},{0,0,1}};

//要素
double element(int LHS_E_or_H,int  LHS_xyz,//未知数の電磁界と方向成分
               int RHS_E_or_H,int  RHS_xyz,//行列のかかる電磁界と方向成分
               int *half,//完全導体境界条件
               int test_x_order,int  test_y_order,
               int basis_x_order,int  basis_y_order,
               int ix ,int iy,
               int x_val,int y_val){
  
  //中心ノード(ix,iy)の周囲のセルの大さを格納
  double 
    Delta_x[2],
    Delta_y[2];
  for(int itex = 0; itex <= 1; itex ++){
    if(ix - 1 + itex < NX[0])
      Delta_x[itex] = DX[0];
    else  if(ix - 1 + itex < NX[1])
      Delta_x[itex] = DX[1];
    else
      Delta_x[itex] = DX[2];
  }
  for(int itey = 0; itey <= 1; itey ++){
    if(iy - 1 + itey < NY[0])
      Delta_y[itey] = DY[0];
    else  if(iy - 1 + itey < NY[1])
      Delta_y[itey] = DY[1];
    else
      Delta_y[itey] = DY[2];
  }
   
  if((LHS_E_or_H == RHS_E_or_H) && (LHS_xyz == RHS_xyz)){//基底関数と基底関数の内積
    int
      ite_start_end[3][2]
      = {
      {0,0},//-1
      {0,1},//0
      {1,1} //1
    };

    /*mu_eps
           itex
          0    1
      1 | 01 | 11 |
 itey   -----n-----
      0 | 00 | 10 |

      n = test_function_center_node
      00,01,10,11 = itex itey
    */
    
    double sum= 0.0;
    for(int itex = ite_start_end[x_val + 1][0]; itex <= ite_start_end[x_val + 1][1]; itex ++){
      for(int itey = ite_start_end[y_val + 1][0]; itey <= ite_start_end[y_val + 1][1]; itey ++){
        
        double eps_mu;
        if(LHS_E_or_H == 0)
          eps_mu
            = dielectrics(
                          position_x( ix - 0.5  + itex*1.0)
                          ,
                          position_y( iy - 0.5  + itey*1.0));
        else
          eps_mu = mu;
        
        sum +=
          eps_mu* //左辺に誘電率,透磁率を乗ずる場合
          inner_product(1
                        ,half[0] ,itex, Delta_x[0], Delta_x[1]
                        ,x_val, test_x_order, basis_x_order)//xに関する内積
          *
          inner_product(1
                        ,half[1] ,itey, Delta_y[0], Delta_y[1]
                        ,y_val, test_y_order, basis_y_order);//yに関する内積
        
      }
    }
    
    return
      sum;
    
  }else if((LHS_E_or_H != RHS_E_or_H) && (LHS_xyz != RHS_xyz)){//基底関数と基底関数の導関数の内積
    double
      eps_mu[2][2];
    double sign = 1.0;
    int minus,plus;
    relative_index(LHS_xyz,&minus,&plus);
    
    if(minus == RHS_xyz) sign = - 1.0;

    
    /*mu_eps
           itex
          0    1
      1 | 01 | 11 |
 itey   -----n-----
      0 | 00 | 10 |

      n = test_function_center_node
      00,01,10,11 = itex itey
    */
    double sum = 0.0;
    
    int
      ite_start_end[3][2]
      = {
      {0,0},//-1
      {0,1},//0
      {1,1} //1
    };
    
    for(int itex = ite_start_end[x_val + 1][0];
        itex <= ite_start_end[x_val + 1][1]; itex ++){
      for(int itey = ite_start_end[y_val + 1][0];
          itey <= ite_start_end[y_val + 1][1]; itey ++){
        
        double eps_mu;
        if(LHS_E_or_H == 0)
          eps_mu
            = dielectrics(
                          position_x( ix - 0.5  + itex*1.0)
                          ,
                          position_y( iy - 0.5  + itey*1.0));
        else
          eps_mu = mu;
        
        sum +=
          inner_product(diff_flag[LHS_xyz][0] || diff_flag[RHS_xyz][0]
                        ,half[0] ,itex, Delta_x[0], Delta_x[1]
                        ,x_val, test_x_order, basis_x_order)//xに関する内積
          *
          inner_product(diff_flag[LHS_xyz][1] || diff_flag[RHS_xyz][1]
                        ,half[1] ,itey, Delta_y[0], Delta_y[1]
                        ,y_val, test_y_order, basis_y_order);//yに関する内積
          // /eps_mu; //右辺に誘電率,透磁率を割る場合。 
          }
    }
  
    return
      sign*DT*sum;
    
  }
  
  std::string str_L_EH ,str_L_xyz ,str_R_EH,str_R_xyz;
  if(LHS_E_or_H == 0)str_L_EH = "電界";
  if(LHS_E_or_H == 1)str_L_EH = "磁界";
  
  if(RHS_E_or_H == 0)str_R_EH = "電界";
  if(RHS_E_or_H == 1)str_R_EH = "磁界";
  
  if(LHS_xyz == 0)str_L_xyz = "x";
  if(LHS_xyz == 1)str_L_xyz = "y";
  if(LHS_xyz == 2)str_L_xyz = "z";

  if(RHS_xyz == 0)str_R_xyz = "x";
  if(RHS_xyz == 1)str_R_xyz = "y";
  if(RHS_xyz == 2)str_R_xyz = "z";
  
  std::cout << str_L_EH  << "の"
            << str_L_xyz << "成分"
            << "を求めるための差分式の右辺に"
            << str_R_EH << "の"
            << str_R_xyz << "成分はありません。"; 

  exit(1);
  
}

void Make_Matrix(Solver **solver
                 ,sparse_matrix ***RHS_Coef_Matrix){
 
  //どの電磁解の成分を求めるかを決定
 
  for(int xyz = 0; xyz <= 2; xyz ++){//0 = Hx, 1 = Hy, 2 = Ez

    int E_or_H = 1;

    if(xyz == 2)E_or_H = 0;
    
    int
      RHS_xyz_minus,
      RHS_xyz_plus; 
    relative_index(xyz,&RHS_xyz_minus,&RHS_xyz_plus);
    int RHS_xyz[3] = {xyz,RHS_xyz_plus,RHS_xyz_minus};
      
    for(int list_index = -1; list_index <= 2; list_index++){
      /*
        Mat[list_index + 1]を生成

        Mat[0] Ez(n+1) = Mat[1] Ez(n) + Mat[2] Hy(n+1/2) + Mat[3] Hx(n+1/2)
       */
        
      int
        target_xyz,
        target_E_or_H = E_or_H;
      if(list_index >= 1)target_E_or_H = !E_or_H;
    
      std::vector <T> list;
    
      if(list_index == -1){
        target_xyz = xyz;
        list.resize((K+1)*10*row_size[xyz]);
      }else{
        
        target_xyz = RHS_xyz[list_index];
        list.resize((K+1)*10*row_size[RHS_xyz[list_index]]);
      }
      if(E_or_H == 1){
        if(list_index > 0 && target_xyz != 2)continue;
      }
      
        
      std::cout<<"Make Coefficient Matrix "<<char_EH[E_or_H]<<char_xyz[xyz]<<list_index + 1<<std::endl;
      //試験関数の中心決定
      for(int x = 0; x <= NX[2]; x++){
        for(int y = 0; y <= NY[2]; y++){
       
            
          int Ncell[2] = {NX[2],NY[2]};
          int position[2] = {x,y};
          int half[2] = {0,0};
          
          //基底関数の片側が完全導体に入っているか判定
          for(int iterator_xyz = 0; iterator_xyz <= 1;  iterator_xyz ++){
            if(0                   == position[iterator_xyz])half[iterator_xyz] = -1;
                
            if(Ncell[iterator_xyz] == position[iterator_xyz])half[iterator_xyz]= 1;
              
          }
              
          //試験関数の導関数の次数 決定
          for(int test_x_order = 0; test_x_order <= K; test_x_order++){
            for(int test_y_order = 0; test_y_order <= K; test_y_order++){
          
              int row = Node( xyz, x, y,  test_x_order, test_y_order);
              if(row == -1)continue;
                    
              //基底関数の中心(試験関数の中心からの距離) 決定
              for(int x_val = -1; x_val <= 1; x_val++){
                for(int y_val = -1; y_val <= 1; y_val++){
   
                  //基底関数の導関数の次数 決定
                  for(int basis_x_order = 0; basis_x_order <= K; basis_x_order++){
                    for(int basis_y_order = 0; basis_y_order <= K; basis_y_order++){
                          
                                
                      int column
                        = Node(target_xyz, x + x_val, y + y_val
                               ,basis_x_order, basis_y_order);
                              
                      if(column != -1){
                           
                        list.push_back(
                                       T(
                                         row,column,
                                         element(
                                                 E_or_H,
                                                 xyz,
                                                 target_E_or_H,
                                                 target_xyz,
                                                 half,
                                                 test_x_order,
                                                 test_y_order,
                                                 basis_x_order,
                                                 basis_y_order,
                                                 x , y,
                                                 x_val, y_val
                                                 )
                                    
                                         )
                                       );
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      
      
      if(list_index == -1){
        sparse_matrix sparse( row_size[xyz],row_size[xyz] );
        sparse.setFromTriplets( list.begin(), list.end() );
        sparse.makeCompressed();
        
        (*solver)[xyz].analyzePattern(sparse);
        (*solver)[xyz].factorize(sparse);
        
      }else{
        (*RHS_Coef_Matrix)[xyz][list_index]
          .resize(row_size[xyz],row_size[target_xyz]);
        
        (*RHS_Coef_Matrix)[xyz][list_index]
          .setFromTriplets(
                           list.begin(),
                           list.end() );
        
        (*RHS_Coef_Matrix)[xyz][list_index]
          .makeCompressed();
      }
    }
  }
 
}










