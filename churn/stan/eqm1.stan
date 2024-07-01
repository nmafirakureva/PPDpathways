data{
  //data points: mean and sd
  real E_m; real E_s;   //remand inflow
  real R_m; real R_s;   //remand stock
  real SI_m; real SI_s; //sentence inflow
  real N_m; real N_s;   //total population
  real FO_a; real FO_b; //fraction open
  real TO_m; real TO_s; //total outflow
  real DR_m; real DR_s; //mean duration at release
  // other priors
  real dR_m; real dR_s;
  real dS_m; real dS_s;
  real dL_m;
  real dom_m;
  real ts_m; real tl_m;
  real ps_a;  real ps_b; 
  real pl_a;  real pl_b; 
  //tolerance for errors
  vector<lower=0>[9] tol;
}
transformed data{
  vector[9] z = rep_vector(0.0,9);
}
parameters{
  //prison populations
  real<lower=0> N;
  real<lower=0> E;
  real<lower=0> R;
  real<lower=0> S;
  real<lower=0> L;
  real<lower=0> Omega;
  // durations/rates
  real<lower=0> dR;
  real<lower=0> dS;
  real<lower=0> dL;
  real<lower=0> dom;
  real<lower=0> ts;
  real<lower=0> tl;
  //proportion
  real<lower=0,upper=1> ps;
  real<lower=0,upper=1> pl;
  // 'data' parameters
  real<lower=0> SI;
  real<lower=0> TO;
  real<lower=0> DR;
  real<lower=0,upper=1> FO;
}
transformed parameters{
  vector[9] eqn; //make vaguely fractional to bring on same scale
  eqn[1] = (E - R/dR)/(1+E);    /* eqm R */
  eqn[2] = (ps*pl*R/dR - tl*L - L/dL) / (1+R); /* eqm L */
  eqn[3] = (ps*(1-pl)*R/dR + tl*L - S/dS - ts*S)/(1+R); /* eqm S */
  eqn[4] = (ts*S - Omega/dom)/(1+Omega);                /* eqm Omega */
  eqn[5] = (ps*R/dR - SI)/(1+R); //proportion sentenced
  eqn[6] = (N - R - S - L - Omega)/(1+N);//total
  eqn[7] = (TO - S/dS - L/dL - Omega/dom)/(1+TO); //total releases NOTE not remand
  /* time served on release from sentence: see notes */
  eqn[8] = ( 1/dR +
             (1/dL)*L/dL +
             ((1/dS + (tl*L/(tl*L+ps*(1-pl)*R/dR))/dL)*S/dS +
              (1/dom + 1/dS + (tl*L/(tl*L+ps*(1-pl)*R/dR))/dL)*Omega/dom)/(1+TO)- DR)/(0.1 + DR);
  eqn[9] = (FO * N - Omega)/(1+N);//fraction in open
}
model{
  //priors for durations
  ps ~ beta(ps_a,ps_b);
  pl ~ beta(pl_a,pl_b);
  dR ~ normal(dR_m,dR_s);
  dS ~ normal(dS_m,dS_s);
  dL ~ exponential(1.0/dL_m);
  dom ~ exponential(1.0/dom_m);
  ts ~ exponential(1.0/ts_m);
  tl ~ exponential(1.0/tl_m);
  // 'data' parameters
  N ~ normal(N_m,N_s);          /* total population */
  E ~ normal(E_m,E_s);          /* remand inflow */
  R ~ normal(R_m,R_s);          /* remand stock */
  SI ~ normal(SI_m,SI_s);       /* sentence inflow */
  TO ~ normal(TO_m,TO_s);       /* total outflow */
  DR ~ normal(DR_m,DR_s);       /* mean duration at release*/
  FO ~ beta(FO_a,FO_b);         /* fraction open */
  //constraint likelihood:
  target += normal_lpdf(eqn | z,tol);
}
