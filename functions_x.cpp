#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double hamiltonian_c(arma::mat P, NumericMatrix Theta, int H_min, int H_max, int L_min, int L_max);
double hamiltonian2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta);
int imax(int a, int b);
int imin(int a, int b);
double dmax(double a, double b);
double dmin(double a, double b);
int getIndex(NumericMatrix potts_info, int h, int l);
int getSum(IntegerMatrix Delta, int k);
IntegerMatrix rectangle2matrix(int H, int L, int H_min, int H_max, int L_min, int L_max, bool rev);
IntegerMatrix potts2ising(IntegerMatrix Delta, int k);
arma::mat potts_c(arma::mat P, NumericMatrix Theta);
arma::mat potts2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta);

// [[Rcpp::export]]
Rcpp::List mcmc_hidden_2(int H, int L, int Q, NumericMatrix potts_info, NumericMatrix Theta_s, NumericMatrix Theta_0_s, IntegerMatrix Delta_s, arma::mat P, double d_s, int iter, int burn, double mu, double mu_0, double sigma, double sigma_0, double dd, double ff, double a_d, double b_d) {
  int i, q, qq, qqq, qqqq, l, h, ll, hh, e, m, count, temp, temp_2;
  int count_2 = 10;
  int M = 1;
  int E = (H*L)*0.05;
  double tau = 0.1;
  double phi = 0.1;
  double hastings = 0;
  double temp_sum, temp_max, d_temp;
  double d = d_s;
  NumericVector accept(4);
  accept(0) = 0;
  accept(1) = 0;
  accept(2) = 0;
  accept(3) = 0;
  NumericVector temp_probs(Q);
  NumericVector temp_count(Q);
  NumericVector d_store(iter);
  NumericMatrix P_1(H, L);
  NumericMatrix P_2(H, L);
  NumericMatrix P_3(H, L);
  NumericMatrix P_4(H, L);
  NumericMatrix theta_store(iter, Q*(Q - 1)/2);
  NumericMatrix theta_0_store(iter, Q*(Q - 1)/2);
  NumericMatrix Theta(Q, Q);
  NumericMatrix Theta_temp(Q, Q);
  NumericMatrix Theta_0(Q, Q);
  arma::mat P_temp(H, L);
  IntegerVector state(Q);
  IntegerMatrix Delta(H, L);
  IntegerMatrix Delta_0(H, L);
  IntegerMatrix Delta_temp(H, L);
  IntegerMatrix Delta_0_temp(H, L);
  IntegerVector Delta_sum(iter);
  NumericMatrix Delta_ppi(H, L);
  
  // Initialization
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
    for (qq = 0; qq < Q; qq++)
    {
      Theta(q, qq) = Theta_s(q, qq);
      Theta_0(q, qq) = Theta_0_s(q, qq);
    }
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      Delta(h, l) = Delta_s(h, l);
      Delta_0(h, l) = 1 - Delta(h, l);
      Delta_ppi(h, l) = 0;
      P_1(h, l) = 0;
      P_2(h, l) = 0;
      P_3(h, l) = 0;
      P_4(h, l) = 0;
      temp_sum = 0;
      for(q = 0; q < Q; q++)
      {
        temp_sum = temp_sum + potts_info(getIndex(potts_info, h + 1, l + 1) - 1, q + 2);
        temp_sum = temp_sum + potts_info(getIndex(potts_info, imax(1, h), imax(1, l)) - 1, q + 2);
        temp_sum = temp_sum + potts_info(getIndex(potts_info, imax(1, h), l + 1) - 1, q + 2);
        temp_sum = temp_sum + potts_info(getIndex(potts_info, h + 1, imax(1, l)) - 1, q + 2);
      }
      if(temp_sum == 0)
      {
        P(h, l) = 0;
        Delta(h, l) = 0;
        Delta_0(h, l) = 1;
      }
    }
  }
  
  // MCMC
  for(i = 0; i < iter; i++)
  {
    // Update d
    d_temp = rgamma(1, d*d/phi, phi/d)(0);
    hastings = (a_d - 1)*(log(d_temp) - log(d)) - b_d*(d_temp - d);
    for(h = 0; h < H; h++)
    {
      for(l = 0; l < L; l++)
      {
        if(P(h, l) != 0)
        {
          if(h != H - 1 && l != L - 1)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P(h, l) - 1) = temp_count(P(h, l) - 1) + 1;
            if(P(h + 1, l) != 0)
            {
              temp_count(P(h + 1, l) - 1) = temp_count(P(h + 1, l) - 1) + 1;
            }
            if(P(h, l + 1) != 0)
            {
              temp_count(P(h, l + 1) - 1) = temp_count(P(h, l + 1) - 1) + 1;
            }
            if(P(h + 1, l + 1) != 0)
            {
              temp_count(P(h + 1, l + 1) - 1) = temp_count(P(h + 1, l + 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d_temp*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings + potts_info(getIndex(potts_info, h + 1, l + 1) - 1, qq + 2)*(d_temp*temp_count(qq) - log(temp_sum));
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings - potts_info(getIndex(potts_info, h + 1, l + 1) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != 0 && l != 0)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P(h, l) - 1) = temp_count(P(h, l) - 1) + 1;
            if(P(h - 1, l) != 0)
            {
              temp_count(P(h - 1, l) - 1) = temp_count(P(h - 1, l) - 1) + 1;
            }
            if(P(h, l - 1) != 0)
            {
              temp_count(P(h, l - 1) - 1) = temp_count(P(h, l - 1) - 1) + 1;
            }
            if(P(h - 1, l - 1) != 0)
            {
              temp_count(P(h - 1, l - 1) - 1) = temp_count(P(h - 1, l - 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d_temp*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings + potts_info(getIndex(potts_info, h, l) - 1, qq + 2)*(d_temp*temp_count(qq) - log(temp_sum));
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings - potts_info(getIndex(potts_info, h, l) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != H - 1 && l != 0)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P(h, l) - 1) = temp_count(P(h, l) - 1) + 1;
            if(P(h + 1, l) != 0)
            {
              temp_count(P(h + 1, l) - 1) = temp_count(P(h + 1, l) - 1) + 1;
            }
            if(P(h, l - 1) != 0)
            {
              temp_count(P(h, l - 1) - 1) = temp_count(P(h, l - 1) - 1) + 1;
            }
            if(P(h + 1, l - 1) != 0)
            {
              temp_count(P(h + 1, l - 1) - 1) = temp_count(P(h + 1, l - 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d_temp*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings + potts_info(getIndex(potts_info, h + 1, l) - 1, qq + 2)*(d_temp*temp_count(qq) - log(temp_sum));
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings - potts_info(getIndex(potts_info, h + 1, l) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != 0 && l != L - 1)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P(h, l) - 1) = temp_count(P(h, l) - 1) + 1;
            if(P(h - 1, l) != 0)
            {
              temp_count(P(h - 1, l) - 1) = temp_count(P(h - 1, l) - 1) + 1;
            }
            if(P(h, l + 1) != 0)
            {
              temp_count(P(h, l + 1) - 1) = temp_count(P(h, l + 1) - 1) + 1;
            }
            if(P(h - 1, l + 1) != 0)
            {
              temp_count(P(h - 1, l + 1) - 1) = temp_count(P(h - 1, l + 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d_temp*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings + potts_info(getIndex(potts_info, h, l + 1) - 1, qq + 2)*(d_temp*temp_count(qq) - log(temp_sum));
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              hastings = hastings - potts_info(getIndex(potts_info, h, l + 1) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
        }
      }
    }
    if (hastings >= log(double(rand()%10001)/10000))
    {
      d = d_temp;
      if(i > burn)
      {
        accept(3)++;
      }
    }
    
    // Update P
    for(e = 0; e < E; e++)
    {
      h = rand()%H;
      l = rand()%L;
      if(P(h, l) != 0)
      {
        for(hh = 0; hh < H; hh++)
        {
          for(ll = 0; ll < L; ll++)
          {
            P_temp(hh, ll) = P(hh, ll);
          }
        }
        for (q = 0; q < Q; q++) 
        {
          P_temp(h, l) = q + 1;
          if(Delta(h, l) == 1)
          {
            temp_probs(q) = hamiltonian2_c(P_temp, Theta, Delta);
          }
          else
          {
            temp_probs(q) = hamiltonian2_c(P_temp, Theta_0, Delta_0);
          }
          if(h != H - 1 && l != L - 1)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P_temp(h, l) - 1) = temp_count(P_temp(h, l) - 1) + 1;
            if(P_temp(h + 1, l) != 0)
            {
              temp_count(P_temp(h + 1, l) - 1) = temp_count(P_temp(h + 1, l) - 1) + 1;
            }
            if(P_temp(h, l + 1) != 0)
            {
              temp_count(P_temp(h, l + 1) - 1) = temp_count(P_temp(h, l + 1) - 1) + 1;
            }
            if(P_temp(h + 1, l + 1) != 0)
            {
              temp_count(P_temp(h + 1, l + 1) - 1) = temp_count(P_temp(h + 1, l + 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              temp_probs(q) = temp_probs(q) + potts_info(getIndex(potts_info, h + 1, l + 1) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != 0 && l != 0)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P_temp(h, l) - 1) = temp_count(P_temp(h, l) - 1) + 1;
            if(P_temp(h - 1, l) != 0)
            {
              temp_count(P_temp(h - 1, l) - 1) = temp_count(P_temp(h - 1, l) - 1) + 1;
            }
            if(P_temp(h, l - 1) != 0)
            {
              temp_count(P_temp(h, l - 1) - 1) = temp_count(P_temp(h, l - 1) - 1) + 1;
            }
            if(P_temp(h - 1, l - 1) != 0)
            {
              temp_count(P_temp(h - 1, l - 1) - 1) = temp_count(P_temp(h - 1, l - 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              temp_probs(q) = temp_probs(q) + potts_info(getIndex(potts_info, h, l) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != H - 1 && l != 0)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P_temp(h, l) - 1) = temp_count(P_temp(h, l) - 1) + 1;
            if(P_temp(h + 1, l) != 0)
            {
              temp_count(P_temp(h + 1, l) - 1) = temp_count(P_temp(h + 1, l) - 1) + 1;
            }
            if(P_temp(h, l - 1) != 0)
            {
              temp_count(P_temp(h, l - 1) - 1) = temp_count(P_temp(h, l - 1) - 1) + 1;
            }
            if(P_temp(h + 1, l - 1) != 0)
            {
              temp_count(P_temp(h + 1, l - 1) - 1) = temp_count(P_temp(h + 1, l - 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              temp_probs(q) = temp_probs(q) + potts_info(getIndex(potts_info, h + 1, l) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
          if(h != 0 && l != L - 1)
          {
            for(qq = 0; qq < Q; qq++)
            {
              temp_count(qq) = 0;
            }
            temp_count(P_temp(h, l) - 1) = temp_count(P_temp(h, l) - 1) + 1;
            if(P_temp(h - 1, l) != 0)
            {
              temp_count(P_temp(h - 1, l) - 1) = temp_count(P_temp(h - 1, l) - 1) + 1;
            }
            if(P_temp(h, l + 1) != 0)
            {
              temp_count(P_temp(h, l + 1) - 1) = temp_count(P_temp(h, l + 1) - 1) + 1;
            }
            if(P_temp(h - 1, l + 1) != 0)
            {
              temp_count(P_temp(h - 1, l + 1) - 1) = temp_count(P_temp(h - 1, l + 1) - 1) + 1;
            }
            temp_sum = 0;
            for(qq = 0; qq < Q; qq++)
            {
              temp_sum = temp_sum + exp(d*temp_count(qq));
            }
            for(qq = 0; qq < Q; qq++)
            {
              temp_probs(q) = temp_probs(q) + potts_info(getIndex(potts_info, h, l + 1) - 1, qq + 2)*(d*temp_count(qq) - log(temp_sum));
            }
          }
        }
        temp_sum = 0;
        temp_max = max(temp_probs);
        for(q = 0; q < Q; q++)
        {
          temp_probs(q) = temp_probs(q) - temp_max;
          temp_probs(q) = exp(temp_probs(q));
          temp_sum = temp_sum + temp_probs(q);
        }
        for(q = 0; q < Q; q++)
        {
          temp_probs(q) = temp_probs(q)/temp_sum;
        }
        temp = RcppArmadillo::sample(state, 1, true, temp_probs)(0) + 1;
        P(h, l) = temp;
      }
    }
    
    // Update Theta
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        for(qqq = 0; qqq < Q; qqq++)
        {
          for (qqqq = 0; qqqq < Q; qqqq++)
          {
            Theta_temp(qqq, qqqq) = Theta(qqq, qqqq);
          }
        }
        Theta_temp(q, qq) = rnorm(1, Theta(q, qq), tau)(0);
        Theta_temp(qq, q) = Theta_temp(q, qq);
        for (h = 0; h < H; h++)
        {
          for (l = 0; l < L; l++)
          {
            P_temp(h, l) = P(h, l);
          }
        }
        for (m = 0; m < M; m++)
        {
          P_temp = potts2_c(P_temp, Theta_temp, Delta);
        }
        hastings = hamiltonian2_c(P_temp, Theta, Delta) - hamiltonian2_c(P, Theta, Delta) + hamiltonian2_c(P, Theta_temp, Delta) - hamiltonian2_c(P_temp, Theta_temp, Delta);
        hastings = hastings - (Theta_temp(q, qq) - mu)*(Theta_temp(q, qq) - mu)/2/sigma/sigma;
        hastings = hastings + (Theta(q, qq) - mu)*(Theta(q, qq) - mu)/2/sigma/sigma;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          Theta(q, qq) = Theta_temp(q, qq);
          Theta(qq, q) = Theta(q, qq);
          if (i > burn) {
            accept(0)++;
          }
        }
      }
    }
    
    // Update Theta_0
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        for(qqq = 0; qqq < Q; qqq++)
        {
          for (qqqq = 0; qqqq < Q; qqqq++)
          {
            Theta_temp(qqq, qqqq) = Theta_0(qqq, qqqq);
          }
        }
        Theta_temp(q, qq) = rnorm(1, Theta_0(q, qq), 0.1)(0);
        Theta_temp(qq, q) = Theta_temp(q, qq);
        for (h = 0; h < H; h++)
        {
          for (l = 0; l < L; l++)
          {
            P_temp(h, l) = P(h, l);
          }
        }
        for (m = 0; m < M; m++)
        {
          P_temp = potts2_c(P_temp, Theta_temp, Delta_0);
        }
        hastings = hamiltonian2_c(P_temp, Theta_0, Delta_0) - hamiltonian2_c(P, Theta_0, Delta_0) + hamiltonian2_c(P, Theta_temp, Delta_0) - hamiltonian2_c(P_temp, Theta_temp, Delta_0);
        hastings = hastings - (Theta_temp(q, qq) - mu_0)*(Theta_temp(q, qq) - mu_0)/2/sigma_0/sigma_0;
        hastings = hastings + (Theta_0(q, qq) - mu_0)*(Theta_0(q, qq) - mu_0)/2/sigma_0/sigma_0;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          Theta_0(q, qq) = Theta_temp(q, qq);
          Theta_0(qq, q) = Theta_0(q, qq);
          if (i > burn) {
            accept(1)++;
          }
        }
      }
    }
    
    // Update Delta
    for(e = 0; e < E; e++)
    {
      h = rand()%H;
      l = rand()%L;
      if (P(h, l) != 0) 
      {
        for(hh = 0; hh < H; hh++)
        {
          for(ll = 0; ll < L; ll++)
          {
            Delta_temp(hh, ll) = Delta(hh, ll);
            Delta_0_temp(hh, ll) = 1 - Delta(hh, ll);
            P_temp(hh, ll) = P(hh, ll);
          }
        }
        Delta_temp(h, l) = 1 - Delta_temp(h, l);
        Delta_0_temp(h, l) = 1 - Delta_0_temp(h, l);
        for (m = 0; m < M; m++)
        {
          P_temp = potts2_c(P_temp, Theta, Delta_temp);
          P_temp = potts2_c(P_temp, Theta_0, Delta_0_temp);
        }
        hastings = hamiltonian2_c(P_temp, Theta, Delta) - hamiltonian2_c(P, Theta, Delta) + hamiltonian2_c(P, Theta, Delta_temp) - hamiltonian2_c(P_temp, Theta, Delta_temp);
        hastings = hastings + hamiltonian2_c(P_temp, Theta_0, Delta_0) - hamiltonian2_c(P, Theta_0, Delta_0) + hamiltonian2_c(P, Theta_0, Delta_0_temp) - hamiltonian2_c(P_temp, Theta_0, Delta_0_temp);
        hh = h;
        ll = l;
        if (Delta(hh, ll) == 0) //Add
        {
          hastings = hastings + (dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))));
        }
        else //Delete
        {
          Delta(hh, ll) = 0;
          hastings = hastings - (dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))));
          Delta(hh, ll) = 1;
        }
        
        if(h != H - 1)
        {
          hh = h + 1;
          ll = l;
          temp = Delta(hh, ll);
          temp_2 = Delta_temp(hh, ll);
          Delta(hh, ll) = 0;
          Delta_temp(hh, ll) = 0;
          hastings = hastings + log(1 + exp(dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))))) - log(1 + exp(dd + ff*(Delta_temp(imin(hh + 1, H - 1), ll) + Delta_temp(imax(hh - 1, 0), ll) + Delta_temp(hh, imin(ll + 1, L - 1)) + Delta_temp(hh, imax(ll - 1, 0)))));
          Delta(hh, ll) = temp;
          Delta_temp(hh, ll) = temp_2;
          if(Delta(hh, ll) == 1)
          {
            if(Delta(h, l) == 0) //Add
            {
              hastings = hastings + ff;
            }
            else //Delete
            {
              hastings = hastings - ff;
            }
          }
        }
        if(h != 0)
        {
          hh = h - 1;
          ll = l;
          temp = Delta(hh, ll);
          temp_2 = Delta_temp(hh, ll);
          Delta(hh, ll) = 0;
          Delta_temp(hh, ll) = 0;
          hastings = hastings + log(1 + exp(dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))))) - log(1 + exp(dd + ff*(Delta_temp(imin(hh + 1, H - 1), ll) + Delta_temp(imax(hh - 1, 0), ll) + Delta_temp(hh, imin(ll + 1, L - 1)) + Delta_temp(hh, imax(ll - 1, 0)))));
          Delta(hh, ll) = temp;
          Delta_temp(hh, ll) = temp_2;
          if(Delta(hh, ll) == 1)
          {
            if(Delta(h, l) == 0) //Add
            {
              hastings = hastings + ff;
            }
            else //Delete
            {
              hastings = hastings - ff;
            }
          }
        }
        if(l != L - 1)
        {
          hh = h;
          ll = l + 1;
          temp = Delta(hh, ll);
          temp_2 = Delta_temp(hh, ll);
          Delta(hh, ll) = 0;
          Delta_temp(hh, ll) = 0;
          hastings = hastings + log(1 + exp(dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))))) - log(1 + exp(dd + ff*(Delta_temp(imin(hh + 1, H - 1), ll) + Delta_temp(imax(hh - 1, 0), ll) + Delta_temp(hh, imin(ll + 1, L - 1)) + Delta_temp(hh, imax(ll - 1, 0)))));
          Delta(hh, ll) = temp;
          Delta_temp(hh, ll) = temp_2;
          if(Delta(hh, ll) == 1)
          {
            if(Delta(h, l) == 0) //Add
            {
              hastings = hastings + ff;
            }
            else //Delete
            {
              hastings = hastings - ff;
            }
          }
        }
        if(l != 0)
        {
          hh = h;
          ll = l - 1;
          temp = Delta(hh, ll);
          temp_2 = Delta_temp(hh, ll);
          Delta(hh, ll) = 0;
          Delta_temp(hh, ll) = 0;
          hastings = hastings + log(1 + exp(dd + ff*(Delta(imin(hh + 1, H - 1), ll) + Delta(imax(hh - 1, 0), ll) + Delta(hh, imin(ll + 1, L - 1)) + Delta(hh, imax(ll - 1, 0))))) - log(1 + exp(dd + ff*(Delta_temp(imin(hh + 1, H - 1), ll) + Delta_temp(imax(hh - 1, 0), ll) + Delta_temp(hh, imin(ll + 1, L - 1)) + Delta_temp(hh, imax(ll - 1, 0)))));
          Delta(hh, ll) = temp;
          Delta_temp(hh, ll) = temp_2;
          if(Delta(hh, ll) == 1)
          {
            if(Delta(h, l) == 0) //Add
            {
              hastings = hastings + ff;
            }
            else //Delete
            {
              hastings = hastings - ff;
            }
          }
        }
        
        if (hastings >= log(double(rand()%10001)/10000))
        {
          Delta(h, l) = 1 - Delta(h, l);
          Delta_0(h, l) = 1 - Delta(h, l);
          if (i > burn) {
            accept(2)++;
          }
        }
      }
    }
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    d_store(i) = d;
    count = 0;
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        theta_store(i, count) = Theta(q, qq);
        theta_0_store(i, count) = Theta_0(q, qq);
        count++;
      }
    }
    temp = 0;
    for(h = 0; h < H; h++)
    {
      for(l = 0; l < L; l++)
      {
        if (i > burn)
        {
          Delta_ppi(h, l) = Delta_ppi(h, l) + Delta(h, l);
          if(P(h, l) == 1)
          {
            P_1(h, l) = P_1(h, l) + 1;
          }
          if(P(h, l) == 2)
          {
            P_2(h, l) = P_2(h, l) + 1;
          }
          if(P(h, l) == 3)
          {
            P_3(h, l) = P_3(h, l) + 1;
          }
          if(P(h, l) == 4)
          {
            P_4(h, l) = P_4(h, l) + 1;
          }
        }
        temp = temp + Delta(h, l);
      }
    }
    Delta_sum(i) = temp;
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      Delta_ppi(h, l) = Delta_ppi(h, l)/(iter - burn);
    }
  }
  accept(0) = accept(0)/(iter - burn)/(Q*(Q - 1)/2);
  accept(1) = accept(1)/(iter - burn)/(Q*(Q - 1)/2);
  accept(2) = accept(2)/(iter - burn)/E;
  accept(3) = accept(3)/(iter - burn);
  return Rcpp::List::create(Rcpp::Named("P_1") = P_1, Rcpp::Named("P_2") = P_2, Rcpp::Named("P_3") = P_3, Rcpp::Named("P_4") = P_4, Rcpp::Named("theta") = theta_store, Rcpp::Named("d") = d_store, Rcpp::Named("theta_0") = theta_0_store, Rcpp::Named("Delta_ppi") = Delta_ppi, Rcpp::Named("Delta_sum") = Delta_sum, Rcpp::Named("accept") = accept);
}

// [[Rcpp::export]]
IntegerMatrix rectangle2matrix(int H, int L, int H_min, int H_max, int L_min, int L_max, bool rev) {
  IntegerMatrix Delta(H, L);
  int h, l;
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if (h >= H_min - 1 & h < H_max & l >= L_min - 1 & l < L_max)
      {
        if(rev) 
        {
          Delta(h, l) = 0;
        }
        else
        {
          Delta(h, l) = 1;
        }
      }
      else
      {
        if(rev) 
        {
          Delta(h, l) = 1;
        }
        else
        {
          Delta(h, l) = 0;
        }
      }
    }
  }
  return Delta;
}

// [[Rcpp::export]]
double hamiltonian_c(arma::mat P, NumericMatrix Theta, int H_min, int H_max, int L_min, int L_max) {
  double hamiltonian = 0;
  //int L = P.n_cols;
  //int H = P.n_rows;
  int h, l;
  for(h = H_min - 1; h < H_max - 1; h++)
  {  
    for(l = L_min - 1; l < L_max - 1; l++)
    {
      if(P(h, l) != 0)
      {
        if(P(h, l) != P(h, l + 1))
        {
          hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
        }
        if(P(h, l) != P(h + 1, l))
        {
          hamiltonian = hamiltonian + Theta(P(h + 1, l) - 1, P(h, l) - 1);
        }
      }
    }
  }
  h = H_max - 1;
  for(l = L_min - 1; l < L_max - 1; l++)
  {
    if(P(h, l) != 0 && P(h, l) != P(h, l + 1))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
    }
  }
  l = L_max - 1;
  for(h = H_min - 1; h < H_max - 1; h++)
  {
    if(P(h, l) != 0 && P(h, l) != P(h + 1, l))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h + 1, l) - 1);
    }
  }
  return (hamiltonian);
}

// [[Rcpp::export]]
double hamiltonian2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta) {
  double hamiltonian = 0;
  int L = P.n_cols;
  int H = P.n_rows;
  int h, l;
  for(h = 0; h < H - 1; h++)
  {  
    for(l = 0; l < L - 1; l++)
    {
      if(P(h, l) != 0 && Delta(h, l) != 0)
      {
        if(Delta(h, l + 1) != 0 && P(h, l) != P(h, l + 1))
        {
          hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
        }
        if(Delta(h + 1, l) != 0 && P(h, l) != P(h + 1, l))
        {
          hamiltonian = hamiltonian + Theta(P(h + 1, l) - 1, P(h, l) - 1);
        }
      }
    }
  }
  h = H - 1;
  for(l = 0; l < L - 1; l++)
  {
    if(P(h, l) != 0 && Delta(h, l) != 0 && Delta(h, l + 1) != 0 && P(h, l) != P(h, l + 1))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
    }
  }
  l = L - 1;
  for(h = 0; h < H - 1; h++)
  {
    if(P(h, l) != 0 && Delta(h, l) != 0 && Delta(h + 1, l) !=0 && P(h, l) != P(h + 1, l))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h + 1, l) - 1);
    }
  }
  return (hamiltonian);
}

// [[Rcpp::export]]
arma::mat potts_c(arma::mat P, NumericMatrix Theta) {
  int L = P.n_cols;
  int H = P.n_rows;
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int h, l, q;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(P(h, l) != 0)
      {
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = 0;
          if(l < L - 1 && (P(h, l + 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l + 1) - 1, q);
          }
          if(l > 0 && (P(h, l - 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l - 1) - 1, q);
          }
          if(h < H - 1 && (P(h + 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h + 1, l) - 1, q);
          }
          if(h > 0 && (P(h - 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h - 1, l) - 1, q);
          }
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat potts2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta) {
  int L = P.n_cols;
  int H = P.n_rows;
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int h, l, q;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(P(h, l) != 0 && Delta(h, l) != 0)
      {
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = 0;
          if(l < L - 1 && Delta(h, l + 1) != 0 && (P(h, l + 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l + 1) - 1, q);
          }
          if(l > 0 && Delta(h, l - 1) != 0 && (P(h, l - 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l - 1) - 1, q);
          }
          if(h < H - 1 && Delta(h + 1, l) != 0 && (P(h + 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h + 1, l) - 1, q);
          }
          if(h > 0 && Delta(h - 1, l) != 0 && (P(h - 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h - 1, l) - 1, q);
          }
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat MRFi(arma::mat P, int H, int L, double e, double f, int iter) {
  int h, l, i;
  double prob;
  for (i = 0; i < iter; i++)
  {
    for (h = 0; h < H; h++)
    {
      for (l = 0; l < L; l++)
      {
        P(h, l) = 0;
        prob = exp(e + f*(P(imin(h + 1, H - 1), l) + P(imax(h - 1, 0), l) + P(h, imin(l + 1, L - 1)) + P(h, imax(l - 1, 0))));
        prob = prob/(1 + prob);
        P(h, l) = rbinom(1, 1, prob)(0);
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat MRFp(arma::mat P, int H, int L, int Q, double f, int iter) {
  int h, l, q, i;
  double temp;
  IntegerVector state(Q);
  NumericVector prob_temp(Q);
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for (i = 0; i < iter; i++)
  {
    for (h = 0; h < H; h++)
    {
      for (l = 0; l < L; l++)
      {
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = 0;
        }
        if (l < L - 1) {
          prob_temp(P(h, l + 1) - 1) = prob_temp(P(h, l + 1) - 1) + f;
        }
        if (l > 0) {
          prob_temp(P(h, l - 1) - 1) = prob_temp(P(h, l - 1) - 1) + f;
        }
        if (h < H - 1) {
          prob_temp(P(h + 1, l) - 1) = prob_temp(P(h + 1, l) - 1) + f;
        }
        if (h > 0) {
          prob_temp(P(h - 1, l) - 1) = prob_temp(P(h - 1, l) - 1) + f;
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
IntegerMatrix potts2ising(IntegerMatrix Delta, int k) {
  int L = Delta.ncol();
  int H = Delta.nrow();
  int l, h;
  IntegerMatrix Delta_k(H, L);
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(Delta(h, l) == k)
      {
        Delta_k(h, l) = 1;
      }
      else
      {
        Delta_k(h, l) = 0;
      }
    }
  }
  return Delta_k;
}

// [[Rcpp::export]]
int getSum(IntegerMatrix Delta, int k) {
  int L = Delta.ncol();
  int H = Delta.nrow();
  int sum = 0;
  int h, l;
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(Delta(h, l) == k)
      {
        sum++;
      }
    }
  }
  return sum;
}

// [[Rcpp::export]]
int getIndex(NumericMatrix potts_info, int h, int l) {
  int index = 0;
  while(h != potts_info(index, 0) || l != potts_info(index, 1))
  {
    index = index + 1;
  }
  return index + 1;
}

// [[Rcpp::export]]
int imax(int a, int b) {
  if(a > b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
int imin(int a, int b) {
  if(a < b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
double dmax(double a, double b) {
  if(a > b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
double dmin(double a, double b) {
  if(a < b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}
