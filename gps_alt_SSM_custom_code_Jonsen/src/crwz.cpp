#include <cmath>
#include <TMB.hpp>

using namespace density; 
using std::sqrt;

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY(Y);
  DATA_VECTOR(tdiff);
  DATA_VECTOR(state0);
  DATA_IVECTOR(isd);                //  indexes observations (1) vs. interpolation points (0)
  DATA_SCALAR(omega);
  DATA_SCALAR(df_z);
  
  // Parameters
  PARAMETER_ARRAY(X);
  PARAMETER_ARRAY(V);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(log_tau);
  PARAMETER(l_rho_p);
  PARAMETER(l_rho_o);
  
  vector<Type> sigma = exp(log_sigma);
  vector<Type> tau = exp(log_tau);
  Type z_err_var = pow(omega,2) + pow(tau(2),2);
  Type z_err_sd = pow(z_err_var, 0.5);
  Type rho_p = Type(2.0) / (Type(1.0) + exp(-l_rho_p)) - Type(1.0);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
//  Type df_z = exp(log_df_z);
  
  // Initial state
  Type jnll = 0.0;
  Type tiny = 1e-5;
  
//  vector<Type> sf = fidx.size();
//  vector<Type> sp = pidx.size();
  
  // Setup object for evaluating multivariate normal likelihood
  matrix<Type> cov(6,6);
  cov.setZero();
  cov(0,0) = tiny;
  cov(1,1) = tiny;
  cov(2,2) = tiny;
  cov(3,3) = sigma(0) * tdiff(0);
  cov(4,4) = sigma(1) * tdiff(0);
  cov(5,5) = sigma(2) * tdiff(0);
  cov(3,4) = pow(sigma(0) * tdiff(0), 0.5) * pow(sigma(1) * tdiff(0), 0.5) * rho_p;
  cov(4,3) = cov(3,4);
  
  // loop over 2 coords and update nll of start location and velocities.
  for(int i = 0; i < Y.rows(); i++) {
    jnll -= dnorm(X(i,0), state0(i), tiny, true);
    jnll -= dnorm(V(i,0), state0(i+3), tiny, true);
  }
  
  // CRW PROCESS MODEL
  vector<Type> x_t(6);
  MVNORM_t<Type> nll_proc(cov); // Multivariate Normal for states
  
  for(int i = 1; i < tdiff.size(); i++) {
    // process cov at time t
    cov.setZero();
    cov(0,0) = tiny;
    cov(1,1) = tiny;
    cov(2,2) = tiny;
    cov(3,3) = sigma(0) * tdiff(i);
    cov(4,4) = sigma(1) * tdiff(i);
    cov(5,5) = sigma(2) * tdiff(i);
    cov(3,4) = pow(sigma(0) * tdiff(i), 0.5) * pow(sigma(1) * tdiff(i), 0.5) * rho_p;
    cov(4,3) = cov(3,4);
    
    // location innovations
    x_t(0) = X(0,i) - (X(0,i-1) + V(0,i) * tdiff(i));
    x_t(1) = X(1,i) - (X(1,i-1) + V(1,i) * tdiff(i));
    x_t(2) = X(2,i) - (X(2,i-1) + V(2,i) * tdiff(i));
    
    // velocity innovations
    x_t(3) = V(0,i) - V(0,i-1); 
    x_t(4) = V(1,i) - V(1,i-1); 
    x_t(5) = V(2,i) - V(2,i-1); 
    
    nll_proc.setSigma(cov);   // set up ith cov matrix
    jnll += nll_proc(x_t);    // Process likelihood
  }
  
  // OBSERVATION MODEL
  // 2 x 2 covariance matrix for observations
  matrix<Type> cov_obs(2, 2);
  MVNORM_t<Type> nll_obs; // Multivariate Normal for observations
  
  
  for(int i = 0; i < tdiff.size(); ++i) {
    if(isd(i) == 1) {
      cov_obs.setZero();
      // GPS observations
      Type x = tau(0); // * K(i,0);
      Type y = tau(1); // * K(i,1);
//      Type z = tau(2);
      cov_obs(0,0) = x * x;
      cov_obs(1,1) = y * y;
      cov_obs(0,1) = x * y * rho_o;
      cov_obs(1,0) = cov_obs(0,1);
 //     cov_obs(2,2) = z * z;
  
      // set up i-th obs & state
      vector<Type> Y_xy(2);
      Y_xy << Y(0,i), Y(1,i);
      vector<Type> X_xy(2);
      X_xy << X(0,i), X(1,i);
  
      nll_obs.setSigma(cov_obs);   // set up i-th obs cov matrix
      jnll += nll_obs(Y_xy - X_xy);   // CRW location innovations
      jnll -= dt((Y(2,i) - X(2,i)) / z_err_sd, df_z, true) - log(z_err_sd);
    }
  }
  

  
  ADREPORT(sigma);
  ADREPORT(rho_p);
  ADREPORT(rho_o);
  ADREPORT(tau);
  REPORT(z_err_sd);
  
  return jnll;
}
