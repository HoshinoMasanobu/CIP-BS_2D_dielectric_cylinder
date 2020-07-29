#include "set_variables.h"
#include "inner_product.h"
#include "make_coefficient_matrix.h"
/*
double inner_product_S(double Delta,int w, int m, int n){
  if(m == 0){
    if(n == 0){
      if(w == 0) return 26.0/35.0 * Delta;
      if(w == -1 ) return 9.0/70.0 * Delta;
      if(w == 1) return 9.0/70.0 * Delta;
    }
    
    if(n == 1){
      if(w == 0) return 0.0;
      if(w == -1 ) return 13.0/420.0 * Delta * Delta;
      if(w == 1)   return -13.0/420.0 * Delta * Delta;
    }
  }

  if(m == 1){
    if(n == 0){
      if(w == 0) return 0.0;
       if(w == -1 ) return -13.0/420.0 * Delta * Delta;
       if(w == 1) return 13.0/420.0 * Delta * Delta;
    }
    
    if(n == 1){
      if(w == 0) return 2.0/105.0 * Delta * Delta * Delta;
      if(w == -1 ) return - 1.0 /140.0 * Delta * Delta * Delta ;
      if(w == 1) return - 1.0 / 140.0 * Delta * Delta * Delta;
    }
  }
  
}
double inner_product_L(double Delta, int w, int m, int n){
  if(m == 0){
    if(n == 0){
      if(w == 0) return 0.0;
      if(w == -1) return - 0.5;
      if(w == 1) return  0.5;
    }
    
    if(n == 1){
      if(w == 0) return 1.0 / 5.0 * Delta;
      if(w == -1 ) return - 1.0 / 10.0 * Delta;
      if(w == 1)   return - 1.0 / 10.0 * Delta;
    }
  }

  if(m == 1){
    if(n == 0){
      if(w == 0) return - 1.0 / 5.0 * Delta;
      if(w == -1 ) return  1.0 / 10.0 * Delta;
      if(w == 1)   return  1.0 / 10.0 * Delta;
    }
    
    if(n == 1){
      if(w == 0) return 0.0;
      if(w == -1 ) return  1.0 /60.0 * Delta * Delta ;
      if(w == 1) return - 1.0 / 60.0 * Delta * Delta;
    }
  }  
}
*/

Eigen::VectorXd *a;
double cf_a(int basis_function_order,int coefficient_order){
  return a[basis_function_order](2*K+1-coefficient_order);
}
double factorial(int n){
  if(n==1||n==0)return 1.0;
  return factorial(n-1)*(double)n;
}
double combination(int n,int m){//nCm
  return factorial(n)/(factorial(n-m)*factorial(m));
}

double inner_product_S(int half,
                       int prev_next,
                       double Delta_prev,double Delta_next,int w, int m, int n){
  
  
  double IP=0.0,tmp;
  
  if(w == 0){
    for(int i = 0; i < K+1; i++){
    
      IP = IP + a[n](i)/((2*K+1-i)+m+1.0) + a[m](i)/((2*K+1-i)+n+1.0);
      for(int h = 0; h < K+1; h++){
    
        IP = IP + a[m](i)*a[n](h)/((4.0*K+2.0-i-h)+1.0); 
      }
    }
    
    double
    tmp_prev = Delta_prev,
    tmp_next = Delta_next;
    
    for(int i = 0; i < m; i++){
      tmp_prev = Delta_prev/((double)m-i)*tmp_prev;
      tmp_next = Delta_next/((double)m-i)*tmp_next;
    }
    for(int i = 0; i < n; i++){
      tmp_prev = Delta_prev/((double)n-i)*tmp_prev;
      tmp_next = Delta_next/((double)n-i)*tmp_next;
    }

    if(half == -1 || prev_next == 1)tmp_prev = 0.0;
    if(half == 1 || prev_next ==0 )tmp_next = 0.0;
    
      return (tmp_next + std::pow(-1.0, m + n)* tmp_prev) * (IP + 1.0/(m+n+1.0));
   
   
  }

  if(half == w){
    //std::cout <<"zero"<<std::endl;
    return 0.0;
  }

  if(prev_next == 0 && w == 1)
    return 0.0;
  if(prev_next == 1 && w == -1)
    return 0.0;
  
  double Delta = Delta_next;
  if(w==-1){
    Delta = Delta_prev;
    tmp = n;
    n = m;
    m = tmp;
  }
  
  tmp = 0.0;
  for(int i = 0; i < K+1; i++){
    for(int l = 0; l <= 2*K+1-i; l++){
      IP = IP 
        + std::pow(-1.0 ,2*K+1-i)*a[n](i)
        * combination(2*K+1-i,l)*std::pow(-1.0, l)
        * 1.0/((2.0*K+1.0-i) + m - l + 1.0);
      for(int j = 0; j < K+1; j++){ 
        IP = IP
          + std::pow(-1.0,2*K+1-i)*a[n](i)
          * combination(2*K+1-i,l)*std::pow(-1.0 ,l)
          * a[m](j)/((4.0*K+2.0-i-j) - l + 1.0); 
      }
    }
  }
  IP=std::pow(-1.0,n)*IP;
  for(int l = 0; l <= n; l++){
    tmp = tmp
      + combination(n,l)* std::pow(-1.0 ,l)
      *1.0/(n - l + m + 1.0);
    for(int j = 0; j < K+1; j++){
      tmp = tmp
        + combination(n,l)* std::pow(-1.0 ,l)
        * a[m](j)/(n - l + (2.0*K + 1.0 - j) + 1.0);
    }
  }
  IP = IP + tmp;
  tmp = Delta;
  for(int i = 0; i < m; i++)tmp = Delta/((double)m-i)*tmp;
  for(int i = 0; i < n; i++)tmp = Delta/((double)n-i)*tmp;
  IP =  IP * tmp;
  return IP;
}
double inner_product_L(int half,
                       int prev_next,
                       double Delta_prev,double Delta_next,
                       int w, int m, int n){

  double IP=0.0,tmp;
  if(w==0){
    
    
    for(int i = 0; i< K+1; i++){
      IP = IP
        + (2.0*K+1.0 - i)*a[n](i)/(2.0*K+1.0 - i + m);//
      for(int j = 0; j<K+1; j++){
        IP = IP + (2*K+1-i)*a[n](i)*a[m](j)/(4.0*K+2.0-i-j);//
      }
    }
    if(n!=0){
      for(int i=0; i<K+1; i++){
        IP = IP
          + n*a[m](i)/(2.0*K+1.0 - i + n);
      }
      IP=IP+(double)n/((double)n+m);
    }
    double
      tmp_prev = 1.0,
      tmp_next = 1.0;
    
    // std::cout<<"check" << tmp_prev <<" " << tmp_next  <<std::endl;
    for(int i = 0; i < m; i++){
      tmp_prev = Delta_prev/((double)m-i)*tmp_prev;
      tmp_next = Delta_next/((double)m-i)*tmp_next;
    }
    for(int i = 0; i < n; i++){
      tmp_prev = Delta_prev/((double)n-i)*tmp_prev;
      tmp_next = Delta_next/((double)n-i)*tmp_next;
    }
    /*
    tmp = 1.0;
    for(int i = 0; i < m; i++)tmp = Delta/((double)m-i)*tmp;
    for(int i = 0; i < n; i++)tmp = Delta/((double)n-i)*tmp;
    */
    if(half == -1 || prev_next == 1)tmp_prev = 0.0;
    if(half == 1 || prev_next == 0)tmp_next = 0.0;
    IP =
      IP *(tmp_next + std::pow(-1.0, m + n + 1.0)*tmp_prev);
    return IP;
  }
  if(half == w) return 0.0;

  
  if(prev_next == 0 && w == 1)
    return 0.0;
  if(prev_next == 1 && w == -1)
    return 0.0;
  

  double Delta = Delta_next;
  if(w==-1){
    Delta = Delta_prev;
    tmp = m;
    m = n;
    n = tmp;
  }
  
  tmp=0.0;
  for(int i = 0; i < K+1; i++){
    for(int l = 0; l <= 2*K+1-i; l++){
      for(int j = 0; j < K+1; j++){
        IP = IP
          + std::pow(-1.0,2*K+1-i)*a[n](i)
          * combination(2*K+1-i,l)*std::pow(-1.0 ,l)
          *(2.0*K+1.0-j)*a[m](j)/((4.0*K+2.0-i-j) - l); 
      }
    }
  }
  for(int l = 0; l <= n; l++){
    for(int j = 0; j < K+1; j++){
      tmp = tmp
        + combination(n,l)* std::pow(-1.0 ,l)
        *(2.0*K+1.0-j)*a[m](j)/(n - l + (2.0*K + 1.0 - j) );
    }
  }
  if(m!=0){
    for(int i = 0; i < K+1; i++){
      for(int l = 0; l <= 2*K+1-i; l++){
        IP = IP 
          + std::pow(-1.0 ,2*K+1-i)*a[n](i)
          * combination(2*K+1-i,l)*std::pow(-1.0, l)
          * (double)m/((2.0*K+1.0-i) + m - l);
      }
    }
    for(int l = 0; l <= n; l++){
      tmp = tmp
        + combination(n,l)* std::pow(-1.0 ,l)
        *(double)m/(n - l + m );
    }
  }
  IP=std::pow(-1.0,n)*IP;
  IP = IP + tmp;
  tmp = 1.0;
  for(int i = 0; i < m; i++)tmp = Delta/((double)m-i)*tmp;
  for(int i = 0; i < n; i++)tmp = Delta/((double)n-i)*tmp;
  IP = -IP*tmp;
  
  if(w==-1)return -IP;
  return IP; 
}

void basis_set(){
  a = new  Eigen::VectorXd [K+1];
  double bb;
  sparse_matrix Psparse(K+1,K+1);
  Eigen::VectorXd b(K+1);
  double* matrixP;
  matrixP = new double [K+1];
  std::vector <T> listP( (K+1)*(K+1) );
  
  for(int i = 0; i < K+1; i++) matrixP[i] = 1.0;
  for(int i = 0; i < K+1; i++){
    for(int h = 0; h < K+1; h++){
     
      listP.push_back( T(i,h, matrixP[h]) );
      matrixP[h]=(2*K+1-i-h) * matrixP[h];  
    }
   
  }

  Eigen::SparseLU < sparse_matrix, Eigen::COLAMDOrdering<int> >
    solverP;
  Psparse.setFromTriplets( listP.begin(), listP.end() );
  Psparse.makeCompressed();
  solverP.analyzePattern(Psparse);
  solverP.factorize(Psparse);
  for(int n = 0; n <= K; n++){
    for(int i = 0; i < K+1; i++)b(i) = 0.0;
    bb=1.0;
    for(int i = 0; i < n+1; i++){
      b(i) = -bb;
      bb=(double)(n-i)*bb;
    }
    a[n] = solverP.solve(b);
  
  }
}
