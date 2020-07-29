#include "set_variables.h"

int relative_index(int i,int *ii/*minus i+1*/,int *iii/*plus i-1*/){
  if(i >= 2)i = -1;
  int j = i + 1;
  *ii = j;
  if(j >= 2)j = -1;
  *iii = j + 1; 
}

void fill_table2(int ******table){
  for(int EH = 0; EH <=1; EH ++){
    for(int xyz = 0; xyz <= 2; xyz++){
      int xyz2,xyz3;
      int order[3],mod[3];
      relative_index(xyz, &xyz2, &xyz3);
      table[xyz][EH][xyz][0][0][0] = EH;
      table[xyz2][EH][xyz][0][0][0] = !EH;
      table[xyz3][EH][xyz][0][0][0] = !EH;
      for(order[xyz] = 0; order[xyz] <= K; order[xyz]++){

        mod[xyz] =  order[xyz]%2;

        for( order[xyz2] = 0; order[xyz2] <= K; order[xyz2] ++){
          mod[xyz2] = order[xyz2]%2;

          for( order[xyz3] = 0; order[xyz3] <= K; order[xyz3] ++){
            mod[xyz3] = order[xyz3]%2;

            for(int ite_mod = 0; ite_mod<=2; ite_mod++){
           
              if(mod[ite_mod] == 0){
                table[ite_mod][EH][xyz][order[0]][order[1]][order[2]]
                  = 
                  table[ite_mod][EH][xyz][0][0][0]; 
              }else if(mod[ite_mod] == 1){
                table[ite_mod][EH][xyz][order[0]][order[1]][order[2]]
                  =
                  !table[ite_mod][EH][xyz][0][0][0];
              }
              
            }
          }
        }
      }
    }
  }
  
}

int fill_table(int max_order,int now_order,int *****table,
               int E_or_H,int x_y_z,int x_diff,int y_diff,int z_diff){
  
  table[E_or_H][x_y_z][x_diff][y_diff][z_diff] = 1;
  int xyz2,xyz3;
  relative_index(x_y_z, &xyz2, &xyz3);
  int index[2][3];
  for(int i = 0;i <= 2; i++){
    index[0][i] = index[1][i] = 0;
  }
  index[0][xyz2] = 1;
  index[1][xyz3] = 1;
  if(max_order == now_order) return 0;
  fill_table( max_order ,now_order+1,  table,
              !E_or_H, xyz3,
              x_diff + index[0][0],
              y_diff + index[0][1],
              z_diff + index[0][2]);
  fill_table( max_order, now_order+1,  table,
              !E_or_H, xyz2
              , x_diff + index[1][0]
              , y_diff + index[1][1]
              , z_diff + index[1][2]);
  return 0;
}
int ******condition_table;
int set_boundary_condition(){
  

  condition_table = new int *****[3];
  
  for(int bound = 0; bound <= 2; bound++){//x = 0 or y = 0
    condition_table[bound] = new int ****[2];
    
    for(int E_or_H = 0; E_or_H <= 1; E_or_H ++){
      condition_table[bound][E_or_H] = new int ***[3];
      
      for(int x_y_z = 0; x_y_z <= 2; x_y_z ++){
        condition_table[bound][E_or_H][x_y_z] = new int **[K + 1];
        
        for(int x_order = 0; x_order <= K; x_order++){
          condition_table[bound][E_or_H][x_y_z][x_order] = new int *[K + 1];
          
          for(int y_order = 0; y_order <= K; y_order ++){
            condition_table[bound][E_or_H][x_y_z][x_order][y_order] = new int [K + 1];

            for(int z_order = 0; z_order <= K; z_order ++){
              condition_table[bound][E_or_H][x_y_z][x_order][y_order][z_order] = 0;
              
            }
          }
        }
      }
    }
  }

 
  
  fill_table2(condition_table);
  return 0;
}

int Boundary_condition(int E_or_H,
                       int x_y_z,
                       int x,int y,
                       int x_order,int y_order){
  if(x == 0 || x == NX[2])
    if(condition_table[0][E_or_H][x_y_z][x_order][y_order][0] == 1)return -1;

  if(y == 0 || y == NY[2])
    if(condition_table[1][E_or_H][x_y_z][x_order][y_order][0] == 1)return -1;


  return 0;
}
int
*****node;
int
*row_size;
void set_node_number(){
  node = new int ****[3];
  
  row_size = new int [3];
  
  for(int x_y_or_z = 0; x_y_or_z <= 2; x_y_or_z ++){
    node[x_y_or_z] = new int ***[NX[2]+1];
      
    for(int x = 0; x <= NX[2]; x ++){
      node[x_y_or_z][x] = new int **[NY[2]+1];
      
      for(int y = 0; y <= NY[2]; y ++){
        node[x_y_or_z][x][y] = new int *[K+1];
            
        for(int x_order = 0; x_order <= K; x_order ++){
          node[x_y_or_z][x][y][x_order] = new int [K+1];
              
          for(int y_order = 0; y_order <= K; y_order ++){
            node[x_y_or_z][x][y][x_order][y_order] = -1;

            
                
          }
        }      
      }
    }
  }

  set_boundary_condition();

  int count = 0 ;
  for(int x_y_or_z = 0; x_y_or_z <= 2; x_y_or_z ++){
  int E_or_H = 1;
  if(x_y_or_z == 2)E_or_H = 0;
  row_size[x_y_or_z] = 0;
  //std::cout << E_or_H << " "<< x_y_or_z <<std::endl;
  for(int x = 0; x <= NX[2]; x ++){
    for(int y = 0; y <= NY[2]; y ++){
     
      for(int x_order = 0; x_order <= K; x_order ++){
        for(int y_order = 0; y_order <= K; y_order ++){
          if(Boundary_condition(E_or_H,
                                x_y_or_z,
                                x,y,
                                x_order,y_order)!=-1){
            //std::cout<<E_or_H<<x_y_or_z<<x<<y<<z<<x_order<<y_order<<z_order<<std::endl;
            node[x_y_or_z][x][y][x_order][y_order]
              = row_size[x_y_or_z];
            row_size[x_y_or_z] ++ ;
                
          }else if(E_or_H == 0 ){
            count++;
          }
        }
      }
    }
  }
 }
  std::cout<<"normal" <<count<<std::endl;
}
int Node(int x_y_or_z, int x, int y,int x_order,int y_order){
  if(x < 0 || x > NX[2]) return -1;
  if(y < 0 || y > NY[2]) return -1;
  if(x_order < 0 || x_order > K)return -1;
  if(y_order < 0 || y_order > K) return -1;

  return node[x_y_or_z][x][y][x_order][y_order];
}

 






































































































































































