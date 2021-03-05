#include "mex.h"
#include "math.h"
#include <armadillo>
#include <time.h>

using namespace arma;
using namespace std;


void matlab2arma_mat(mat& A, const mxArray *mxdata){
  // delete [] A.mem; // don't do this!
  access::rw(A.mem)=mxGetPr(mxdata);
  access::rw(A.n_rows)=mxGetM(mxdata); // transposed!
  access::rw(A.n_cols)=mxGetN(mxdata);
  access::rw(A.n_elem)=A.n_rows*A.n_cols;
};

void matlab2arma_vec(vec& A, const mxArray *mxdata){
  // delete [] A.mem; // don't do this!
  access::rw(A.mem)=mxGetPr(mxdata);
  access::rw(A.n_rows)=mxGetM(mxdata); // transposed!
  access::rw(A.n_cols)=mxGetN(mxdata);
  access::rw(A.n_elem)=A.n_rows*A.n_cols;
};


void free_mat(mat& A, const double *ptr){
  access::rw(A.mem)=ptr;
  access::rw(A.n_rows)=1; // transposed!
  access::rw(A.n_cols)=1;
  access::rw(A.n_elem)=1;
};

void free_vec(vec& A, const double *ptr){
  access::rw(A.mem)=ptr;
  access::rw(A.n_rows)=1; // transposed!
  access::rw(A.n_cols)=1;
  access::rw(A.n_elem)=1;
};



// Hier beginnt der essentielle Code

void timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    printf("%s",asctime( localtime(&ltime) ) );
}

double delta_infinity( double x, double y){
   double result; 
   if (std::abs(x-y)<pow(10,-15)){
        result=0;

   }
   else{
       result=std::max(x,y);
   }
  return(result);
}
mat construct_ultra_cost_mat_one(const mat& ux, const mat& uy, const mat& coupling){
  int n = ux.n_cols;
  int m = uy.n_cols;
  mat cost_mat;
  cost_mat.zeros(n,m);

  /*std::cout <<"ux:"<<ux<<std::endl;
  std::cout <<"uy:"<<uy<<std::endl;
  std::cout <<"coupling:"<<coupling<<std::endl;
  std::cout <<"p:"<<p<<std::endl;*/




for (size_t i=0; i<n; i++){
    for (size_t j=0; j<m; j++){
        double tmp=0;
        for (size_t k=0; k<n; k++){
            for (size_t l=0; l<m; l++){
                tmp=tmp+delta_infinity(ux.at(i,k),uy.at(j,l))*coupling.at(k,l);
               // std::cout <<"delta_infinity(ux[i,k],uy[j,l])^p:"<<pow(delta_infinity(ux[i,k],uy[j,l]),p)<<std::endl;
            }
        }
        cost_mat.at(i,j)=2*tmp;
        /*std::cout <<"cost_mat(i,j):"<<cost_mat[i,j]<<std::endl;
        std::cout <<"(i):"<<i<<std::endl;
        std::cout <<"(j):"<<j<<std::endl;*/


    }
}
          
  return(cost_mat);
}

// Hier endet der essentielle Code. 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //timestamp();
  if (nrhs != 3)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");
  
  mat ux(1,1);
  const double* ux_mem=access::rw(ux.mem);
  matlab2arma_mat(ux,prhs[0]);
  
  mat uy(1,1);
  const double* uy_mem=access::rw(uy.mem);
  matlab2arma_mat(uy,prhs[1]);
  
  mat coupling(1,1);
  const double* coupling_mem=access::rw(coupling.mem);
  matlab2arma_mat(coupling,prhs[2]);
    
  mxArray *output = mxCreateDoubleMatrix(ux.n_rows, uy.n_cols, mxREAL);
  mxDouble *out_mem = mxGetPr(output);
  
  //timestamp();
  
  mat cost_mat(1,1);  
  cost_mat = construct_ultra_cost_mat_one(ux, uy, coupling);
  //timestamp();

  const double* cost_mat_mem=access::rw(cost_mat.mem);

  //std::cout <<"cost_mat:"<<cost_mat<<std::endl;
  std::memcpy(out_mem, cost_mat_mem, sizeof(double)*cost_mat.n_elem);  
  plhs[0]=output;
  
  free_mat(ux,ux_mem);
  free_mat(uy,uy_mem);
  free_mat(coupling,coupling_mem);
  free_mat(cost_mat,cost_mat_mem);
  //timestamp();

  return;
  } 
