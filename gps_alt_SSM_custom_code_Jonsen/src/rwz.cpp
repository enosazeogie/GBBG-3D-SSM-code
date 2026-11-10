#include <cmath>
#include <TMB.hpp>

using namespace density; 
using std::sqrt;

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0) {
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY(Y);
  DATA_VECTOR(tdiff);
  DATA_VECTOR(state0);
  DATA_IVECTOR(isd);                //  indexes observations (1) vs. interpolation points (0)
  DATA_SCALAR(df_z);
  DATA_SCALAR(omega);
  
  // Parameters
  PARAMETER_ARRAY(X);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(log_tau);
//  PARAMETER(log_df_z);
  
  vector<Type> sigma = exp(log_sigma);
  vector<Type> tau = exp(log_tau);
  Type z_err_var = pow(omega,2) + pow(tau(2),2);
  Type z_err_sd = pow(z_err_var, 0.5);
//  Type df_z = exp(log_df_z);
  
  // Initial state
  Type jnll = 0.0;
  Type tiny = 1e-5;
  
  // update nll of start location
  jnll -= dnorm(X(0,0), state0(0), tiny, true);
  jnll -= dnorm(X(0,1), state0(1), tiny, true);
  jnll -= dnorm(X(0,2), state0(2), tiny, true);
  
  // RW PROCESS MODEL
  for(int i = 1; i < tdiff.size(); i++) {
    // location innovations
    jnll -= dnorm(X(i,0), X(i-1,0), sigma(0) * tdiff(i), true);
    jnll -= dnorm(X(i,1), X(i-1,1), sigma(1) * tdiff(i), true);
    jnll -= dnorm(X(i,2), X(i-1,2), sigma(2) * tdiff(i), true);
  }
  
  // OBSERVATION MODEL
  for(int i = 0; i < tdiff.size(); ++i) {
    if(isd(i) == 1) {
    // GPS observations
    jnll -= dnorm(Y(i,0), X(i,0), tau(0), true);
    jnll -= dnorm(Y(i,1), X(i,1), tau(1), true);
    jnll -= dt((Y(i,2) - X(i,2)) / z_err_sd, df_z, true) - log(z_err_sd);
//    jnll -= dnorm(Y(i,2), X(i,2), omega, true);
    } 
  }
  
  ADREPORT(sigma);
  ADREPORT(tau);
  REPORT(z_err_sd);
//  ADREPORT(df_z);
  
  return jnll;
}
