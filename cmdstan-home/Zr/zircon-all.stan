//T LaBone
//September 2, 2021
//612 model

//Start FUNCTIONS
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
  
int count_sub(int[] sub,int index) {
  int count;
  count = 0;
  for(i in 1:num_elements(sub)) {
    if(sub[i] == index) {
        count = count + 1;
    }
  }
  return(count);
}

int count_bld(int[] sub,vector dT,int index) {
  int count;
  count = 0;
  for(i in 1:num_elements(sub)) {
    if(sub[i] == index) {
      if(dT[i] == 0.0) {
        count = count + 1;
      }
    }
  }
  return(count);
}

int count_urn(int[] sub,vector dT,int index) {
  int count;
  count = 0;
  for(i in 1:num_elements(sub)) {
    if(sub[i] == index) {
      if(dT[i] > 0.0) {
        count = count + 1;
      }
    }
  }
  return(count);
}

int[] inds(int[] sub,int index) {
  int res[count_sub(sub,index)];
  int ci;
  ci = 1;
  for(i in 1:num_elements(sub)) {
    if(sub[i] == index) {
      res[ci] = i;
      ci = ci + 1;
    }
  }
  return(res);
}

int[] indu(int[] sub,vector dT,int index) {
  int res[count_urn(sub,dT,index)];
  int ci;
  ci = 1;
  for(i in 1:num_elements(sub)) {
    if(sub[i]==index) {
      if(dT[i] > 0.0) {
        res[ci] = i;
        ci = ci + 1;
      }
    }
  }
  return(res);
}

int[] indb(int[] sub,vector dT,int index) {
  int res[count_bld(sub,dT,index)];
  int ci;
  ci = 1;
  for(i in 1:num_elements(sub)) {
    if(sub[i]==index) {
      if(dT[i] == 0.0) {
        res[ci] = i;
        ci = ci + 1;
      }
    }
  }
  return(res);
}
  
} //End FUNCTIONS 
 
// Start DATA
data {
  int ntot;
  int N;
  int n[N];
  int sub[ntot];
  vector[ntot] T;
  vector[ntot] dT;
  vector[ntot] M;
  vector[ntot] u;
  int nt;
  int nc;
  int H[nt,2];
  vector[N] mu_v;
  vector[N] sigma_v;
  vector[nt] mu_theta;
  matrix[nt,nt] sigma_theta;
  real mu_tau;
  real sigma_tau;
  real Eta;
  real df;
} // End DATA

//PARAMETERS
parameters {
  vector[N] V;
  matrix[nt,N] Kt;
  corr_matrix[nt] Omega;
  real<lower=0> tau;
  vector[nt] theta;
}
  
//TRANSFORMED PARAMETERS
transformed parameters{
  cov_matrix[nt] Pi; 
  vector[nt] Tau;
  vector[nt] kt;
  real v;
  vector[ntot] m;

  Tau = rep_vector(tau,nt);
  Pi = quad_form_diag(Omega,Tau);
  
  for(i in 1:N) {
    kt = exp(Kt[1:nt,i]);   
    v = exp(V[i]);
    
    m[indb(sub,dT,i)] = (
      q_c(1, T[indb(sub,dT,i)], kt, H, nc, nt) +
      q_c(2, T[indb(sub,dT,i)], kt, H, nc, nt)) / v;
      
    m[indu(sub,dT,i)] = q_c(15, T[ indu(sub,dT,i)], kt, H, nc, nt) -
      q_c(15, T[ indu(sub,dT,i)] - dT[ indu(sub,dT,i)], kt, H, nc, nt);
      
    }
}

//MODEL
model {
  Omega ~ lkj_corr(Eta);
  theta ~ multi_normal(mu_theta,sigma_theta);
  tau ~ normal(mu_tau,sigma_tau);
  
  for(i in 1:N) {
    Kt[1:nt,i] ~ multi_normal(theta,Pi);
    V[i] ~  normal(mu_v[i],sigma_v[i]);
    M[inds(sub,i)] ~ student_t(df,m[inds(sub,i)],u[inds(sub,i)]);
  }
  
}

//Start GENERATED QUANTITIES
generated quantities{
  real M_hat[ntot]; 
  real log_lik[ntot];
  for(i in 1:N) {
    M_hat[inds(sub,i)] =  student_t_rng(df,m[inds(sub,i)],u[inds(sub,i)]);
    for(j in inds(sub,i)) {
      log_lik[j] = student_t_lpdf(M[j] | df,m[j],u[j]);
    }
  }
}//End Generated QUANTITIES




