#ifndef FOURIE_TRANSFORM
#define FOURIE_TRANSFORM

void set_fourie_trans_parametar();
void update_FT(Eigen::VectorXd **EM_density,int nt);
void output_result_FourieTrans();


void set_fourie_trans_parametar_circle();
void update_FT_circle(Eigen::VectorXd **EM_density,int nt);
void output_result_FourieTrans_circle();

#endif

