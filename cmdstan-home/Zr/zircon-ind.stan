//T LaBone
//August 14, 2021
//612 model

  //FUNCTIONS
  functions {
  
  //Calculate total removal rate constants in matrix K
  matrix trrc( matrix k, real lambda) {
    real Ksum;
    matrix[rows(k),rows(k)] K;
    K = k;
    for(i in 1:rows(K)) {
      Ksum = 0;
      for(j in 1:rows(K)){
        Ksum = Ksum + k[i,j];
      }
      K[i,i] = -Ksum - lambda;
    }
    return K;
  }
  
  //calculate compartment content using the matrix exponential method
  vector q_me(int c,vector x, vector q0, matrix k) {
    matrix[rows(k),rows(k)] a;
    vector[rows(x)] q; 
    q = rep_vector(0,rows(x));
    for (i in 1:rows(x)) {
      a = matrix_exp(k*x[i]);
      for (j in 1:rows(k)) {  
        q[i] = q[i] + a[j,c]*q0[j];
      }
    }
    return q;
  }
  
  //calculate content at times x and rate constants kt
  vector q_c( int comp, vector x, vector kt, int[,] H, int nc, int nt ) {
    vector[nc] q0;
    matrix[nc,nc] k;
    vector[rows(x)] content;
    q0 = rep_vector(0,nc);
    q0[1] = 1.0;
    k = rep_matrix(0,nc,nc);
    for(i in 1:nt) {
      k[H[i,1],H[i,2]] = kt[i];
    }
    k = trrc(k,0);
    content = q_me(comp,x,q0,k);
    return content ;
  }
  
  } 
  
  // DATA
  data {
    int<lower=1> np;
    int<lower=1> nu;
    int<lower=1> n;
    vector[n] M;
    vector[n] u;
    vector[n] ru;
    vector[n] T;
    vector[n] dT;
    int<lower=1> nt;
    int<lower=1> nc;
    int H[nt,2];
    vector[nt] mu_theta;
    matrix[nt,nt] sigma_theta;
    real mu_tau;
    real sigma_tau;
    real mu_beta;
    real sigma_beta;
    real mu_v;
    real sigma_v;
    real Eta;
    real df;
  }

  //PARAMETERS
  parameters {
    real beta;
    corr_matrix[nt] Omega;
    vector[nt] theta;
    real<lower=0> tau;
    real v;
    vector[nt] Kt;
  }
  
  //TRANSFORMED PARAMETERS
  transformed parameters{
    vector[np] m_p;
    vector[nu] m_u;
    vector[n] m;
    vector[nt] kt;
    real Beta;
    real V;
    cov_matrix[nt] Pi; 
    vector[nt] Tau;
    vector[n] mu;
    vector[n] sigma;
    
    Tau = rep_vector(tau,nt);
    Pi = quad_form_diag(Omega,Tau);
    
    kt = exp(Kt);
    Beta = exp(beta);
    V = exp(v);
    m_p = ( q_c(1,T[1:np],kt,H,nc,nt) + q_c(2,T[1:np],kt,H,nc,nt) ) / V;
    m_u = q_c(15,T[(np+1):n],kt,H,nc,nt) - q_c(15,T[(np+1):n]-dT[(np+1):n],kt,H,nc,nt);
    m = append_row(m_p,m_u);
    sigma = u; 
    mu = Beta * m;
  }

  //MODEL
  model {
    Omega ~ lkj_corr(Eta);
    theta ~ multi_normal(mu_theta,sigma_theta);
    tau ~ normal(mu_tau,sigma_tau);
    v ~  normal(mu_v,sigma_v);
  
    Kt ~ multi_normal(theta,Pi);
    beta ~ normal(mu_beta,sigma_beta);
    M ~ student_t(df,mu,sigma);
  }

  //GENERATED QUANTITIES
  generated quantities{
    vector[n] M_hat;
    for(i in 1:n) {
      M_hat[i] = student_t_rng(df,Beta*m[i],Beta*m[i]*ru[i]); 
    }
  }

