functions{
  real binormal_cdf(real z1, real z2, real rho) {
  if (z1 != 0 || z2 != 0) {
    real denom = fabs(rho) < 1.0 ? sqrt((1 + rho) * (1 - rho)) : not_a_number();
    real a1 = (z2 / z1 - rho) / denom;
    real a2 = (z1 / z2 - rho) / denom;
    real product = z1 * z2;
    real delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2);
  }
  return log((0.25 + asin(rho) / (2 * pi())));
}
  //real icar_normal_lpdf(vector phi, int N, array[ ] int node1, array[ ] int node2){
   // return -0.5*dot_self(phi[node1] - phi[node2]);
 // }
  real my_custom_likelihood_lpdf(vector y, matrix x, matrix z, real alpha, real beta10, real beta20,  int N,  vector beta1, vector beta2, vector phi1, vector phi2, vector phi3, int[] K, vector beth, real betc){
    int t = size(y);
    // int f = t/2;
    int f = t/2;
    vector[f] y1;
    vector[f] y2;
    for (i in 1:f){
      y1[i] = y[i];
      int ff = f+i;
      y2[i] = y[ff]; 
    }
    real loglike = 0;
    real loglike1 = 0;
    real loglike2 = 0;
    real loglike3 = 0;
    real loglike4 = 0;
    //real rho
      for (i in 1:N) {
         int k = K[i];
        real rho  = tanh(x[i]*beth+betc + phi3[K[i]]);
        //for (k in 1: 37){
        loglike1 = loglike+((1-y1[i])*(1-y2[i]))* (binormal_cdf(Phi_approx((-x[i]*beta1)-beta10 - phi1[k]), Phi_approx((-z[i]*beta2)-beta20 - phi2[k]), -rho ));
        loglike2 = loglike+ ((y1[i])*(1-y2[i]))* (binormal_cdf(1, Phi_approx((-z[i]*beta2)-beta20 - phi2[k]), rho) - binormal_cdf(Phi_approx((-x[i]*beta1)-beta10 - phi1[k]), Phi_approx((-z[i]*beta2)-beta20 - phi2[k]), -rho));
        loglike3 = loglike+((1-y1[i])*(y2[i]))* (binormal_cdf(Phi_approx((-x[i]*beta1)-beta10 - alpha - phi1[k]), 1, -rho) - binormal_cdf(Phi_approx((-x[i]*beta1)-beta10 - phi1[k]), Phi_approx((-z[i]*beta2)-beta20 - phi2[k]), rho));
        loglike4 = loglike+((y1[i])*(y2[i]))* (1-binormal_cdf(Phi_approx((-x[i]*beta1)-beta10-alpha - phi1[k]), 1, rho) - binormal_cdf(1,Phi_approx((-z[i]*beta2-beta20) - phi2[k]), rho) + binormal_cdf(Phi_approx((-x[i]*beta1)-beta10 - alpha - phi1[k]), Phi_approx((-z[i]*beta2)-beta20 - phi2[k]),rho));
          //}
        }
      loglike = loglike1 + loglike2 + loglike3 + loglike4;
      return loglike;
    }
}
data {
  int <lower=0> N; // number of data items
  //int t = 2*N;
  vector[2 * N] y;       // Observed data (dependent variable)
  matrix[N,22] x;
  matrix[N,22] z;
  //int<lower=0, upper=k> N_edges;
  int K[N];
  vector[37] mu; 
  //matrix[37, 37] scale_matrix_inv;
  //real<lower=0> df;
  matrix[37, 37] D_W_inv; 
  //matrix[37,37] D_W_inv;
  //array[N_edges] int<lower=1, upper=N>node1; //node1[i] adjacent to node2[i]
  //array[N_edges] int<lower=1, upper=N>node2; // and node1[i] < node2[i]
  //int<lower=0,upper=1> y[N];
}
parameters { 
  vector[22] beta1;
  vector[22] beta2;
  real beta10;
  real beta20;
  vector[37] phi1; 
  vector[37] phi2;
  vector[37] phi3;
  vector[22] beth;
  //vector[37] mu; 
  //real Sigma; 
  real betc;
  //real<lower=0> sigma1;
  //real<lower=0> sigma2;
  real alpha;
  //real rho;
  //matrix[2,2] Sigma;
  //real nu;
  }
//transformed parameters{
//  vector[N] rho;
//  vector[N] h;
//  for (i in 1:N){
    //h[i] = x[i]*beth+betc + phi[K[i]];
//    rho[i] = tanh(x[i]*beth+betc + phi[K[i]]);
 // }
//}
model{
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1); 
  beth ~ normal(0,1);
  betc ~ normal(0,1);
  //rho ~ normal(0,1000);
  //h ~ normal(0,1000);
  alpha ~ normal(0,1);
  beta10 ~ normal(0,1);
  beta20 ~ normal(0,1);
  //Sigma ~ inv_gamma(0.4,0.14);
  //mu ~ normal(0,100);
  //Sigma ~ inv_wishart(df, scale_matrix_inv);
  phi1 ~ multi_normal(mu, 10*D_W_inv);
  phi2 ~ multi_normal(mu, 10*D_W_inv);
  phi3 ~ multi_normal(mu, 10*D_W_inv);
  y ~ my_custom_likelihood(x, z, alpha, beta10, beta20, N, beta1, beta2, phi1, phi2, phi3, K,beth,betc); 
}
