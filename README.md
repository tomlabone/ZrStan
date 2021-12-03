# ZrStan
Stan code to perfom hierarchical Bayesian analysis of zirconium bioassay data
===========

Introduction
------------

All of the data preparation and calculations in my dissertation were performed with R scripts and cmdstan code (Carpenter et al. 2017) run from bash scripts. While completely reproducible, this body of work is complex and difficult to communicate to others concisely. Therefore, to provide additional details on the calculations presented in Chapters *3* and *4* of the dissertation, here I generate multiple sets of simulated bioassay data for *16* hypothetical subjects. The simulated data from some subset of the *16* subjects (e.g., *5* subjects to speed up the calculations) are then used as the input data for the population model calculations discussed in Chapter *3*. The canned priors calculated from the simulated population data are then used to evaluate the simulated bioassay data from a couple individuals not used in development of the canned priors, which is akin to the calculations performed in Chapter 4. Finally, a description of the codes and procedures are provided. The goal is to allow interested individuals to download and run the codes, which will facilitate them reproducing and expanding on what I have done.

Simulated Data
--------------

Simulated bioassay data were calculated starting with the sex, weight, *t*, *Δt*, *u*, *n*, *n*<sub>u</sub>, and *n*<sub>p</sub> from the *16* subjects as defined in Chapter *3*. For each subject *100* prior predictive distributions were generated using the canned priors for *θ* and *τ* from Chapter 4 and the blood plasma volume posterior distributions from Chapter 3. After the simulated blood plasma data were generated the blood plasma volume for each subject was calculated using the methods described in Appendix B of the disseration. For example, some bioassay data for Subject Zr\_01 from run *50* of the simulation are shown in the table below and spaghetti plots of the data from all *16* subjects for run *50* are shown in the plots. Data with *Δt=0* are blood plasma samples and those with *Δt>0* are urine samples The units of the blood plasma bioassay are the fraction of a unit quantity of injected zirconium per liter of blood plasma and the units of the incremental urine bioassay are fraction of a unit quantity of injected zirconium present in the urine sample collected from *t-Δt* to *t*. Weight has units of kg, vol liters, *t* days, and *Δt* days. These data are in the file *ZrSimData-612.csv*.

| Run  |  Subject |  Sex |   Weight |    Vol  |         *t*    |            *Δt*  |         M   |               u |
|:----:| :--------:| :---:| :--------:| :---:|     :------------:|    :---:  |:--------:| :-----------------:|
|50    | Zr\_01   |    m  |      81 |    4.956 |  0.0034722 |     0   |      0.1981    |  0.018209 |
| 50   |   Zr\_01 |    m  |      81 |    4.956 |  0.013889   |     0   |      0.18535  |     0.0092612 |
| 50   |  Zr\_01  |    m  |      81 |    4.956 |  0.020139   |     0   |      0.19931  |     0.0085196 |
| 50   |   Zr\_01 |    m  |      81 |    4.956 |  0.030556   |     0   |      0.19551  |     0.0083008 |




![Simulated zirconium urine and blood plasma data for *16* subjects from run *50* of the simulation.[\[fig:Simulated-zirconium-urine\]]{#fig:Simulated-zirconium-urine label="fig:Simulated-zirconium-urine"}](UrineSpagSim.png "fig:")
![Simulated zirconium urine and blood plasma data for *16* subjects from run *50* of the simulation.[\[fig:Simulated-zirconium-urine\]]{#fig:Simulated-zirconium-urine label="fig:Simulated-zirconium-urine"}](PlasmaSpagSim.png "fig:")



The R script *ZrDataSim.R* reads in the data from a specified run of the simulated data, the *50th* run in this example, and creates a data file *Zr612-50.data.R* and an init file *Zr612-50.init.R* that are readable by the cmdstan code *zircon-all.stan*.

Computing Environment
---------------------

The calculations in this chapter were performed on two computers running Linux Mint 20.2 having 128Gb of RAM, one of which had an i7-7820x cpu and the other an i9-10900X cpu. The statistical computing software R[^1] was run in the integrated development environment RStudio[^2] that was used to prepare datasets for and analyze the output from the program cmdstan[^3] that did the MCMC calculations. Bash scripts were used to execute the cmdstan code as a batch file running in the background of the Linux operating system. Each MCMC chain was run on its own core. I tried disabling hyper threading to speed up the calculations but found that the calculations ran faster with hyper threading on.

Although Stan was used in my dissertation, the MCMC software package NIMBLE (Valpine et al., 2017) was also considered. NIMBLE is very promising but it did not offer methods for solving the system of differential equations, i.e., the matrix exponential function was not available. NIMBLE did offer eigenvalue decomposition, which can also be used to solve the systems of differential equations, but it tended to return imaginary numbers for rate matrices (which is not physically possible) and was therefore impractical to use. NIMBLE is designed to be extensible, and my future plans include implementing the matrix exponential and eigenvalue decomposition method so that its performance can be compared to Stan's. Stan did not offer an eigenvalue decomposition method for non-symmetric matrices like those encountered with biokinetic models, and I also plan on extending Stan to handle these calculations.

Overview
--------

The following steps give instructions on how to assemble simulated bioassay datasets for populations and individuals, analyze them with Stan using the cmdstan interface, and evaluate the posterior distributions with R code.

1.  Run *ZrDataSim.R* to generate cmdstan input data files for the *16* individuals from run *50* of the simulated dataset and cmdstan input data files for the first *5* of these individuals who are designated as the study population in this example.

2.  Run cmdstan program *zircon-all.stan* from the bash script *Zr-all-50.sh*.

3.  Run *Zranal-Sim.R* to analyze Stan output and create canned prior dataset.

4.  Run the first part of *Zrind-Sim.R* to prepare the cmdstan data file for subject *10* from run *50* of the simulated data and to perform a standard regression of the reference bioassay functions on the simulated bioassay data for this subject.

5.  Run cmdstan program *zircon-in.stan* from the bash script *Zr\_10.sh*.

6.  Run the second half of *Zrind-Sim.R* to analyze the Stan output.

The three R scripts are in the main level of the GitHub directory and are not discussed any further because they are only used to prepare and analyze cmdstan output data. The sub directories of GitHub are

-   cmdstan-home: bash scripts (and when you run it, cmdstan)

    -   Zr: cmdstan code

-   data: input data

-   diag: MCMC diagnostic plots

-   plots: output plots

-   subjects: datasets for individuals from run 50 of the simulated data.

The cmdstan code and bash scripts are discussed in more detail in the following sections.

Bash script for Population Model
--------------------------------

The bash script *Zr-all-50.sh* shown below runs the cmdstan code with *3000* samples for the warm up followed by *10000* samples from the posterior distribution. The biokinetic model parameters, prior distributions, and bioassay data are read from *Zr612-50.data.R* and the starting values for the parameters from *Zr612-50.init.R*. Note that if a starting value is not defined for a given parameter Stan uses a randomly drawn value. The posterior samples for each parameter are written to *Zr612-50.csv*. This script assumes that it and the folder Zr reside in the default cmdstan-home folder as described in the documentation for cmdstan [@stan_team_cmdstan_2021].

    #!/bin/bash
    Zr/zircon-all \
       sample num_warmup=3000 num_samples=10000 \
       data file=Zr/Zr612-50.data.R \
       init=Zr/Zr612-50.init.R \
       output file=Zr/Zr612-50.csv

The command used to submit the bash script as a background job is

    nohup ./Zr-all-50.sh > Zr-50.out &

where the warnings, errors, and progress information are written to *Zr-50.out*.

Stan Code for Population Model
------------------------------

Cmdstan code resembles a mix of C++ and R, and consists of the following sections:

-   Functions - user defined functions.

-   Data - anything that has a fixed value and is not subject to stochastic variation.

-   Parameters - anything that is assigned a probability distribution in the Model section.

-   Transformed Parameters - anything that takes a parameter and performs a calculation with it before drawing a sample from the posterior distribution. Transformed parameter affect the sampling of the posterior.

-   Model - the declaration of the probability distributions of the parameters.

-   Generated Quantities - anything that takes a parameter and performs a calculation with it after drawing a sample from the posterior distribution. Generated quantities do not affect the sampling of the posterior.

The cmdstan code *zircon-all.stan* for the zirconium population model discussed in this section must be compiled as described in the cmdstan documentation before use. Note that all data structures used in cmdstan code must be declared before use.

### Functions

The code in the Function section is broken down into two parts, the first giving functions used to solve the system of ordinary differential equations and the second to track and assign data and parameters to each of the 16 individuals. The key function in the first part is *q\_c,* which calculates reference bioassay functions given the compartment number *comp*, times *x*, rate matrix *kt*, adjacency matrix *H*, number of compartments *nc*, and, number of non-zero rate constants *nt* that are read from *Zr612-50.data.R.* The adjacency matrix concisely represents how the nodes in a digraph like the one used to represent the ICRP 134 zirconium model are connected (Birchall, 1989).

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

The bioassay results for all subjects are arranged in one data vector and there is another vector that contains a number from *1* to *16* that indicates the subject to whom each bioassay result belongs. The following functions parse the data vector into blood plasma and urine data for a given subject. The data were arranged this way because each subject can have a different number of bioassay results and Stan does not support ragged matrices.

    //gives the number of bioassay results for a given subject   
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
    
    //gives the number of blood plasma results for a given subject
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
    
    //gives the number of urine results for a given subject
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
    
    //returns the indices for all results for a given subject
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
    
    //returns the indices for urine results for a given subject
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
    
    //returns the indices for urine results for a given subject
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
    }}

### Data

The biokinetic model and bioassay data are read into the cmdstan code from *Zr612-50.data.R.*

    data {
      int ntot;                  //total number of bioassay results
      int N;                     //number of individuals
      int n[N];                  //number of results for each individual
      int sub[ntot];             //the subject number for each bioassay result
      vector[ntot] T;            //the time of each bioassay result
      vector[ntot] dT;           //the delta t for each bioassay result
      vector[ntot] M;            //the bioassay results
      vector[ntot] u;            //the measurement uncertainties
      int nt;                    //the number of transfer rate constants
      int nc;                    //the number of compratments
      int H[nt,2];               //the adjacency matrix
      vector[N] mu_v;            //the log mean of the plasma volume
      vector[N] sigma_v;         //the log sd of the plasma volume
      vector[nt] mu_theta;       //log means of theta
      matrix[nt,nt] sigma_theta; //log covariance matrix of theta
      real mu_tau;               //log mean of tau
      real sigma_tau;            //log sd of tau
      real Eta;                  //Eta for LKJ distribution
      real df;                   //degrees of freedom for Student t
    }

### Parameters

Parameters have probability density models that are defined in the Model section.

    parameters {
      vector[N] V;
      matrix[nt,N] Kt;
      corr_matrix[nt] Omega;
      real tau;
      vector[nt] theta;
    }

### Transformed Parameters

In the Transformed Parameter section values of the parameters drawn in the Model section are manipulated to give variables that are a function of one or more parameters. The data functions are first used in this section. For example, `T[indb(sub,dT,i)]` retrieves the times of all blood plasma measurements for the ith subject.

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
          q_c(15, T[indu(sub,dT,i)] - dT[ indu(sub,dT,i)], kt, H, nc, nt);
        }
    }

### Models

The Model section is where probability distributions are assigned.

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

### Generated Quantities.

In the Generated Quantities section the current sample of the posterior distributions of the parameters are used to generate the predicted bioassay results for the posterior predictive plot and the log likelihood in case one wants to compare different models.

    generated quantities{
      real M_hat[ntot]; 
      real log_lik[ntot];
      for(i in 1:N) {
        M_hat[inds(sub,i)] =  student_t_rng(df,m[inds(sub,i)],u[inds(sub,i)]);
        for(j in inds(sub,i)) {
          log_lik[j] = student_t_lpdf(M[j] | df,m[j],u[j]);
        }
      }
    }

Bash script for Individual Model
--------------------------------

The bash script *Zr\_10.sh* shown below contains the *sample* command that is used to calculate the posterior distributions of the parameters in the individual model (analogous to the bash script for the population model) and the *optimize* command that is used to generate MAP estimates of the parameters. One nice thing about Stan is that the same code is used for both calculations, i.e., we don't have to write any new code to get the MAP estimates.

    Zr/zircon-ind \
       sample num_warmup=4000 num_samples=15000 \
       data file=Zr/Zr_10.data.R \
       init=Zr/Zr_10.init.R \
       output file=Zr/Zr_10.csv
       
    Zr/zircon-ind \
      optimize \
      data file=Zr/Zr_10.data.R \
      init=Zr/Zr_10.init.R \
      output file=Zr/Zr_10opt.csv

Stan Code for Individual Model
------------------------------

The Stan code for the individual model is somewhat simpler than the code for the population model because it deals with a single individual, which means that the *6* functions in used to parse the bioassay data for multiple subjects in the population model are not needed and there are no loops in the Model and Transformed Parameters sections of the code.

### Functions

The functions to solve the system of ordinary differential equations are the same as those in the population model.

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
    }} 

### Data

The data are for one subject, but otherwise are essentially the same as in the Data section of the population code.

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

### Parameters

The intake *β* for an individual subject is a parameter to estimate unlike in the population model where it is a known constant.

    parameters {
      real beta;
      corr_matrix[nt] Omega;
      vector[nt] theta;
      real<lower=0> tau;
      real v;
      vector[nt] Kt;
    }

### Transformed Parameters

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

### Model

The all the subjects in the population model had known quantities of zirconium injected that were normalized to *β=1* in order to facilitate comparisons of bioassay data from subject to subject. When the data from an individual are modeled the quantity injected *\beta* is an unknown parameter to be estimated. Having the intake and blood plasma volume as parameters to estimate is a significant source of uncertainty in the final parameter estimates for a out-of-sample subject.

    model {
      Omega ~ lkj_corr(Eta);
      theta ~ multi_normal(mu_theta,sigma_theta);
      tau ~ normal(mu_tau,sigma_tau);
      v ~  normal(mu_v,sigma_v);
      Kt ~ multi_normal(theta,Pi);
      beta ~ normal(mu_beta,sigma_beta);
      M ~ student_t(df,mu,sigma);
    }

### Generated Quantities

The main difference between the Generated Quantities sections for the population versus the individual model is that *σ* for the t distribution is estimated from the product of the predicted bioassay measurement `Beta*m[i]` times the relative uncertainty in the observed bioassay measurement `ru[i] = u[i]/M[i]`.

    generated quantities{
      vector[n] M_hat;
      for(i in 1:n) {
        M_hat[i] = student_t_rng(df,Beta*m[i],Beta*m[i]*ru[i]); 
      }
    }

Procedure
---------

To run the population model:

1.  Compile *zircon-all.stan*, noting that it is assumed to be in the cmdstan-home/Zr directory.

2.  Place *Zr612-50.data.R* and *Zr612-50.init.R.* in the Zr directory.

3.  Place *Zr-all-50.sh* in the cmdstan-home directory, set its permissions, and submit it as a background job using

    a.  chmod 700 Zr-all-50.sh

    b.  nohup ./Zr-all-50.sh \> Zr-all-50.out &

4.  Note that the output goes to the text file *Zr-all-50.out* and it takes about *4* days for the example to run using the first *5* subjects as the population. In comparison, analyzing all *16* subjects takes about *25* days.

5.  Run the R code *Zranal-Sim.R* to perform diagnostics and generate canned priors.

To run the individual model:

1.  Compile *zircon-ind.stan*, noting that it is assumed to be in the cmdstan-home/Zr directory.

2.  Run the first part of the R code Zrind-Sim.R, which reads in the canned priors from *Zrcan-50.R* and the bioassay data for the 10th subject from *Zr\_10.R*.

3.  Place *Zr\_10.data.R* and *Zr\_10.init.R.* in the Zr directory.

4.  Place *Zr-ind-10.sh* in the cmdstan-home directory, set its permissions, and submit it as a background job using

    a.  chmod 700 Zr-ind-10.sh

    b.  nohup ./Zr-ind-10.sh \> Zr-ind-10.out &

5.  Note that the output goes to the text file *Zr-ind-10.out* and it takes about *3* hours for the example to run.

6.  Run the second part of the R code *Zrind-Sim.R* to perform diagnostics and generate output.

[^1]: Version 4.1 available from https://cran.r-project.org/

[^2]: Version 1.4.1717 available from https://www.rstudio.com/products/rstudio/

[^3]: Version V2.27 available from https://mc-stan.org/users/interfaces/cmdstan

Bibliography
------
Birchall, A. and A. C. James (June 1989). “A microcomputer algorithm for solving first-order compartmental models involving recycling.” In: Health physics 56.6, pp. 857–868. issn: 0017-9078. doi: 10.1097/00004032-198906000-00003.

Carpenter, Bob et al. (Jan. 11, 2017). “Stan: A Probabilistic Programming Language”. In: Journal of Statistical Software 76.1. Number: 1, pp. 1–32. issn: 1548-7660. doi: 10.18637/jss.v076.i01.

Valpine, Perry de et al. (Apr. 3, 2017). “Programming With Models:Writing Statistical Algorithms for General Model Structures With NIMBLE”. In: Journal of Computational and Graphical Statistics 26.2, pp. 403–413. issn:1061-8600. doi:10.1080/10618600.2016.1172487.
