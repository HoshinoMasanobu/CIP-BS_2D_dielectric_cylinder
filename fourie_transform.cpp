#include "set_variables.h"
#include "tfsf.h"
#include "interpolation.h"
#include <complex>
#include <cmath>
using Complex = std::complex<double>;

double FT_DX,FT_DY;
int FT_NX = 100,FT_NY = 100;
class vector2d{
public:
  int ix,iy;
};
vector2d *iterator_FT;
Complex *F_Ez, *F_Ez_broadband;

vector2d *iterator_FT_circle;
Complex *F_Ez_circle, *F_Ez_broadband_circle;

int max_iterator,max_iterator_circle;
void set_fourie_trans_parametar(){
  iterator_FT = new vector2d[2*(FT_NX + FT_NY)];
  F_Ez = new Complex[2*(FT_NX + FT_NY)];
  F_Ez_broadband = new Complex[NT+1];
  
  int
    update_xy[4][2] = {{1,0},{0,1},{-1,0},{0,-1}};
  
  FT_DX = (FT_x_end - FT_x_start)/(double)FT_NY;
  FT_DY = (FT_y_end - FT_y_start)/(double)FT_NX;
  std::cout<<FT_DX<<" "<<FT_DY<<std::endl;
  int
    FT_corner[4][2] = {{0,0},{FT_NX,0},{FT_NX,FT_NY},{0,FT_NY}};
  
  int iez = 0;
  
  for(int corner = 0; corner <= 3; corner ++){
       
    int
      next_corner = (corner + 1)%4;
    
    
    int ix = FT_corner[corner][0] ;
    do{
      int iy = FT_corner[corner][1];
      do{
          
        iterator_FT[iez].ix = ix;
        iterator_FT[iez].iy = iy;
        F_Ez[iez] = 0.0;
        iez++;
      
        iy += update_xy[corner][1];
      }while( iy != FT_corner[next_corner][1] );
    
      ix += update_xy[corner][0];
    }while( ix != FT_corner[next_corner][0] );
  }
  max_iterator = iez;
  std::cout << max_iterator<<std::endl;
}

void update_FT( Eigen::VectorXd **EM_density, int nt){
  Complex zj(0.0,1.0);
  
  for(int i = 0; i < max_iterator; i++  ){
    
    double
      x = FT_DX*iterator_FT[i].ix + FT_x_start + DX[0]*pml_L,
      y = FT_DY*iterator_FT[i].iy + FT_y_start + DY[0]*pml_L;
    
    F_Ez[i] +=
      DT*EM_interpolated(EM_density,2, x, y)
      *std::exp(- zj * omega * DT * (double)nt);
    
  }
  /*
  double
    x = FT_DX*iterator_FT[200].ix + FT_x_start + DX[0]*pml_L,
    y = FT_DY*iterator_FT[200].iy + FT_y_start + DY[0]*pml_L;
  
  for(int i = 0; i <= NT; i++){
    double freq = 2.0*M_PI * i/(double)NT;

    F_Ez_broadband[i] +=
      DT*EM_interpolated(EM_density,2, x, y)
      *std::exp(- zj * freq * (double)nt);
  }
  */
}

void output_result_FourieTrans(){
  std::string filename = "FourieTrans_Ez_NX"
    +std::to_string(NX[0] - pml_L)+"_"+std::to_string(NX[1] - NX[0])+".dat";
  std::ofstream ofs(filename);
  for(int i = 0; i < max_iterator; i++  ){
    
    double
      x = FT_DX*iterator_FT[i].ix + FT_x_start,
      y = FT_DY*iterator_FT[i].iy + FT_y_start;

    ofs << x << " " << y << " "<< FT_DX*i << " "
        << std::abs(F_Ez[i]) << " "
        << std::atan2(F_Ez[i].imag(),F_Ez[i].real()) << std::endl;
  }
  ofs.close();

  std::string filename2 = "FT_broadband.dat";
  std::ofstream ofs2(filename2);
 
  for(int i = 0; i <= NT; i++){
    double freq = 2.0*M_PI/(double)NT * i/DT;

    ofs2 << freq << " " << std::abs(F_Ez_broadband[i]) <<" "
         << std::atan2(F_Ez_broadband[i].imag(),F_Ez_broadband[i].real()) << std::endl;
   
  }
  ofs2.close();
  std::cout<<omega<<std::endl;
}

void set_fourie_trans_parametar_circle(){

  F_Ez_circle = new Complex[200];
  F_Ez_broadband_circle = new Complex[NT+1];
  /*
  int
    update_xy[4][2] = {{1,0},{0,1},{-1,0},{0,-1}};
  
  FT_DX = (FT_x_end - FT_x_start)/(double)FT_NY;
  FT_DY = (FT_y_end - FT_y_start)/(double)FT_NX;
  */ 
  //std::cout<<FT_DX<<" "<<FT_DY<std::endl;
  int iez = 0;

  max_iterator_circle = 200;
}

void update_FT_circle( Eigen::VectorXd **EM_density, int nt){
  Complex zj(0.0,1.0);
  double
    LX = RX[0] + RX[1] + RX[2],
    LY = RY[0] + RY[1] + RY[2];
  double Radius = std::sqrt(2.0)*(1.0/2.0 + 1.0/2.0);
  for(int i = 0; i < max_iterator_circle; i++  ){
    
    double
      x = Radius*std::cos(M_PI + phi_i + 2.0*M_PI*i/(double)max_iterator_circle)  + LX/2.0,
      y = Radius*std::sin(M_PI + phi_i + 2.0*M_PI*i/(double)max_iterator_circle)  + LY/2.0;
    
    F_Ez_circle[i] +=
      DT*EM_interpolated(EM_density,2, x, y)
      *std::exp(- zj * omega * DT * (double)nt);
    
  }
  
  double
    x = Radius*std::cos(phi_i )  + LX/2.0,
    y = Radius*std::sin(phi_i )  + LY/2.0;
  
  for(int i = 0; i <= NT; i++){
    double freq = 2.0*M_PI * i/(double)NT;

    F_Ez_broadband_circle[i] +=
      DT*EM_interpolated(EM_density,2, x, y)
      *std::exp(- zj * freq * (double)nt);
  }
}


void output_result_FourieTrans_circle(){
 double
    LX = RX[0] + RX[1] + RX[2],
    LY = RY[0] + RY[1] + RY[2];
  double Radius = std::sqrt(2.0)*(1.0/2.0 + 1.0/2.0);
  std::string filename = "FourieTrans_Ez_NX"
    +std::to_string(NX[0] - pml_L)+"_"+std::to_string(NX[1] - NX[0])
    +"Phi"+std::to_string(phi_i)+".dat";
  std::ofstream ofs(filename);
  for(int i = 0; i < max_iterator_circle; i++  ){
    
    double
      x = Radius*std::cos(M_PI + phi_i + 2.0*M_PI*i/(double)max_iterator_circle)  + LX/2.0,
      y = Radius*std::sin(M_PI + phi_i + 2.0*M_PI*i/(double)max_iterator_circle)  + LY/2.0;
    
    ofs << x << " " << y << " "<< 2.0*M_PI*i/(double)max_iterator_circle << " "
        << std::abs(F_Ez_circle[i]) << " "
        << std::atan2(F_Ez_circle[i].imag(),F_Ez_circle[i].real()) << std::endl;
  }
  ofs.close();

  std::string filename2 = "FT_broadband_circle.dat";
  std::ofstream ofs2(filename2);
  
  for(int i = 0; i <= NT; i++){
    double freq = 2.0*M_PI/(double)NT * i/DT;
    
    ofs2 << freq << " " << std::abs(F_Ez_broadband_circle[i]) <<" "
         << std::atan2(F_Ez_broadband_circle[i].imag(),F_Ez_broadband_circle[i].real()) << std::endl;
    
  }
  ofs2.close();
  std::cout<<omega<<std::endl;
}


































