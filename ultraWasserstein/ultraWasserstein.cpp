#include "mex.h"
#include "math.h"
#include <armadillo>

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

double one_path( arma::uvec index_set, const vec& mu,const vec& nu, const mat& dist_mat,const double& p,int& position,  std::vector<uvec>& global_indexsets,std::vector<vec>& global_distances, std::vector<int>& global_realsubset,std::vector<double>& global_height){
  mat dist_mat_changing = dist_mat.submat(index_set,index_set);
  vec dist_vec= global_distances[position];
  //std::cout << "position:"<< position <<std::endl;
  //std::cout << "dist_vec:"<< dist_vec <<std::endl;
  double value =0;
  if(index_set.size()==1){
    //std::cout << "1.1" <<std::endl;
    value = value + abs(sum(mu.elem(index_set)-nu.elem(index_set)))*pow(global_height[position]/2,p);
    //std::cout << value <<std::endl;
  }
  else{
    for (size_t i=0; i<(dist_vec.size()-1); ++i){
      uvec set=index_set.elem(find(dist_mat_changing.col(0)<dist_vec[i]));
      uvec set_complement= index_set.elem(find(dist_mat_changing.col(0)>=dist_vec[i]));
      //  std::cout << "i:"<<i <<std::endl;
      //  std::cout << "SET:"<<set <<std::endl;
      //  std::cout << "SET_COM:"<<set_complement <<std::endl;
      
      position = position +1;
      
      global_indexsets[position]= set_complement;
      global_height[position]=dist_vec[i];
      
      mat tmp_one = dist_mat.submat(set_complement,set_complement);
      vec tmp_two =tmp_one.col(0);
      arma::vec dist_vec_set_complement = reverse(unique(tmp_two));
      
      global_distances[position]= dist_vec_set_complement;
      
      if(dist_vec_set_complement[0]<dist_vec[i]){
        global_realsubset[position]=1;
      }
      
      value = value + abs(sum(mu.elem(set)-nu.elem(set)))*(pow(dist_vec[i]/2,p)-pow(dist_vec[i+1]/2,p));
      /*std::cout << "value:"<<value <<std::endl;
       std::cout << "summandmu:"<<abs(sum(mu.elem(set)-nu.elem(set))) <<std::endl;
       std::cout << "summanddist:"<<(pow(dist_vec[i]/2.,p)-pow(dist_vec[i+1]/2.,p)) <<std::endl;
       std::cout << "summanddist1:"<<(pow(dist_vec[i]/2.,p)) <<std::endl;
       std::cout << "summanddist2:"<<(pow(dist_vec[i+1]/2.,p)) <<std::endl;
       std::cout << "dist_vec:"<<dist_vec <<std::endl;
       std::cout << "set:"<<set <<std::endl;*/
      
      
      dist_mat_changing = dist_mat.submat(set,set);
      index_set = set;
    }
  }
  return(value);
}

double one_path_real_subset( arma::uvec index_set, const vec& mu,const vec& nu, const mat& dist_mat,const double& p,int& position,  std::vector<uvec>& global_indexsets,std::vector<vec>& global_distances, std::vector<int>& global_realsubset,std::vector<double>& global_height){
  mat dist_mat_changing = dist_mat.submat(index_set,index_set);
  vec dist_vec= global_distances[position];
  
  double value =0;
  uvec set = index_set;
  
  if(index_set.size()==1){
    value = value +abs(sum(mu.elem(index_set)-nu.elem(index_set)))*pow(global_height[position]/2,p);
    /*std::cout << "value:"<<value <<std::endl;
     std::cout << "summandmu:"<<abs(sum(mu.elem(set)-nu.elem(set))) <<std::endl;
     std::cout << "summanddist:"<<pow(global_height[position]/2,p)<<std::endl;
     std::cout << "dist_vec:"<<dist_vec <<std::endl;
     std::cout << "indexset:"<<index_set <<std::endl;*/
  }
  else{
    for (size_t i=0; i<(dist_vec.size()); ++i){
      
      value = value + abs(sum(mu.elem(set)-nu.elem(set)))*(pow(global_height[position]/2,p)-pow(dist_vec[i]/2,p));
      /*std::cout << "value:"<<value <<std::endl;
       std::cout << "summandmu:"<<abs(sum(mu.elem(set)-nu.elem(set))) <<std::endl;
       std::cout << "summanddist:"<<(pow(dist_vec[i]/2.,p)-pow(dist_vec[i+1]/2.,p)) <<std::endl;
       std::cout << "summanddist1:"<<(pow(global_height[position]/2,p)) <<std::endl;
       std::cout << "summanddist2:"<<(pow(dist_vec[i]/2,p)) <<std::endl;
       std::cout << "dist_vec:"<<dist_vec <<std::endl;
       std::cout << "set:"<<set <<std::endl;*/
      set=index_set.elem(find(dist_mat_changing.col(0)<dist_vec[i]));
      
      if(dist_mat_changing.n_cols>1){
        position = position +1;
        uvec set_complement= index_set.elem(find(dist_mat_changing.col(0)>=dist_vec[i]));
        
        global_indexsets[position]= set_complement;
        global_height[position]=dist_vec[i];
        
        mat tmp_one = dist_mat.submat(set_complement,set_complement);
        vec tmp_two =tmp_one.col(0);
        arma::vec dist_vec_set_complement = reverse(unique(tmp_two));
        
        global_distances[position]=dist_vec_set_complement;
        
        /*std::cout << "set_complement:"<<set_complement <<std::endl;
         std::cout << "dist_vec_set_complement:"<<dist_vec_set_complement <<std::endl;
         double global_heighttmp = global_height[position-1];
         std::cout << "global_height:"<<global_heighttmp <<std::endl;*/
        
        if(dist_vec_set_complement[0]<global_height[position-1]){
          global_realsubset[position]= 1;
        }
        dist_mat_changing = dist_mat.submat(set,set);
        index_set = set;
      }
    }
  }
  return(value);
}
// [[Rcpp::export]]
double UltrametricWasserstein(const vec& mu, const vec& nu, const mat& dist_mat, double& p){
  int n = mu.size();
  std::vector<vec> global_distances(n);
  std::vector<uvec> global_indexsets(n);
  
  
  std::vector<int> global_realsubset(n, 0);
  std::vector<double> global_height(n, 0);
  double ultraWasser = 0;
  int position = 0;
  
  global_indexsets[position]= linspace<uvec>(0, n-1,n);
  global_distances[position]= reverse(unique(dist_mat.col(0)));
  
  while(position >=0){
    int position_old = position;
    if(global_realsubset[position]==0){
      ultraWasser=ultraWasser+one_path(global_indexsets[position],mu, nu, dist_mat,p,position,global_indexsets,global_distances,global_realsubset,global_height);
      /*std::cout << "ultraWasser:"<<ultraWasser <<std::endl;
       std::cout << "global_realsubset:" <<std::endl;
       for (auto i = global_realsubset.begin(); i !=global_realsubset.end(); ++i){
       std::cout << *i << ' ';
       }
       std::cout << "\n" <<std::endl;*/
      
    }
    else{
      ultraWasser = ultraWasser+one_path_real_subset(global_indexsets[position],mu, nu, dist_mat,p,position,global_indexsets,global_distances,global_realsubset,global_height);
      /*std::cout << "ultraWasser:"<<ultraWasser <<std::endl;
       std::cout << "global_realsubset:" <<std::endl;
       for (auto i = global_realsubset.begin(); i !=global_realsubset.end(); ++i){
       std::cout << *i << ' ';
       }
       std::cout << "" <<std::endl;*/
    }
    //std::cout <<"position:"<<position<<std::endl;
    //uvec tmp_1 = global_indexsets[1];
    //vec tmp_2 = global_distances[1];
    //std::cout <<"global_indexsets:"<<tmp_1<<std::endl;
    //std::cout <<"global_distances:"<<tmp_2<<std::endl;
    
    global_indexsets.erase(global_indexsets.begin()+position_old);
    global_distances.erase(global_distances.begin()+position_old);
    
    
    global_height.erase(global_height.begin()+position_old);
    global_realsubset.erase(global_realsubset.begin()+position_old);
    position = position -1;
  }
  ultraWasser = pow(2,(p-1))*ultraWasser;
  return(ultraWasser);
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 4)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");
  
  vec mu(1);
  const double* mu_mem=access::rw(mu.mem);
  matlab2arma_vec(mu,prhs[0]);
  
  vec nu(1);
  const double* nu_mem=access::rw(nu.mem);
  matlab2arma_vec(nu,prhs[1]);
  
  mat dist_mat(1,1);
  const double* dist_mat_mem=access::rw(dist_mat.mem);
  matlab2arma_mat(dist_mat,prhs[2]);
  
  double p = mxGetScalar(prhs[3]);
  
  double result = UltrametricWasserstein(mu, nu, dist_mat, p);
  
  plhs[0] = mxCreateDoubleScalar((double)result);

  free_vec(mu,mu_mem);
  free_vec(nu,nu_mem);
  free_mat(dist_mat,dist_mat_mem);
  return;
  } 
