# include <Rcpp.h>

using namespace Rcpp ;

// [[Rcpp::export()]]

NumericVector get_trigonometric_summand_3_parameter(NumericVector x, int J_1, int J_2, int J_3) {
  
  int i,j_1,j_2,j_3,I;
  NumericVector T(1+2*J_1*J_2*J_3);
  double c, s, c_1, s_1, c_2, c_3, m_1, m_2, m_3;
  
  for (i = 0; i < x.size(); i++){
    T[0] += 1;
    c_1 = cos(x[i]);
    c_2 = cos(2*J_1*x[i]);
    c_3 = cos(2*J_2*J_1*x[i]);
    m_1 = c_1*c_1;
    c = c_1;
    s_1 = sqrt(1-m_1)*((x[i] < 0) ? -1 : (x[i] > 0));
    s = s_1;
    m_2 = 1;
    m_3 = 1;
    I = 1;
    
    for (j_3 = 0; j_3 < J_3-1; j_3++){
      for (j_2 = 0; j_2 < J_2; j_2++){
        for (j_1 = 0; j_1 < J_1; j_1++){
          T[I] += c;
          c *= m_1;
          I += 1;
          T[I] += s;
          s *= m_1;
          I += 1;
        }
        m_2 *= c_2;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      m_3 *= c_3;
      m_2 = m_3;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_2 = 0; j_2 < J_2-1; j_2++){
      for (j_1 = 0; j_1 < J_1; j_1++){
        T[I] += c;
        c *= m_1;
        I += 1;
        T[I] += s;
        s *= m_1;
        I += 1;
      }
      m_2 *= c_2;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_1 = 0; j_1 < J_1-1; j_1++){
      T[I] += c;
      c *= m_1;
      I += 1;
      T[I] += s;
      s *= m_1;
      I += 1;
    }
    T[I] += c;
    I += 1;
    T[I] += s;
  }
  
  return T;
}

// [[Rcpp::export()]]

NumericVector get_trigonometric_summand_4_parameter(NumericVector x, int J_1, int J_2, int J_3, int J_4) {
  
  int i,j_1,j_2,j_3,j_4,I;
  NumericVector T(1+2*J_1*J_2*J_3*J_4);
  double c, s, c_1, s_1, c_2, c_3, c_4, m_1, m_2, m_3, m_4;
  
  for (i = 0; i < x.size(); i++){
    T[0] += 1;
    c_1 = cos(x[i]);
    c_2 = cos(2*J_1*x[i]);
    c_3 = cos(2*J_2*J_1*x[i]);
    c_4 = cos(2*J_3*J_2*J_1*x[i]);
    m_1 = c_1*c_1;
    c = c_1;
    s_1 = sqrt(1-m_1)*((x[i] < 0) ? -1 : (x[i] > 0));
    s = s_1;
    m_2 = 1;
    m_3 = 1;
    m_4 = 1;
    I = 1;
    
    for (j_4 = 0; j_4 < J_4-1; j_4++){
      for (j_3 = 0; j_3 < J_3; j_3++){
        for (j_2 = 0; j_2 < J_2; j_2++){
          for (j_1 = 0; j_1 < J_1; j_1++){
            T[I] += c;
            c *= m_1;
            I += 1;
            T[I] += s;
            s *= m_1;
            I += 1;
          }
          m_2 *= c_2;
          c = c_1*m_2;
          s = s_1*m_2;
        }
        m_3 *= c_3;
        m_2 = m_3;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      m_4 *= c_4;
      m_3 = m_4;
      m_2 = m_3;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_3 = 0; j_3 < J_3-1; j_3++){
      for (j_2 = 0; j_2 < J_2; j_2++){
        for (j_1 = 0; j_1 < J_1; j_1++){
          T[I] += c;
          c *= m_1;
          I += 1;
          T[I] += s;
          s *= m_1;
          I += 1;
        }
        m_2 *= c_2;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      m_3 *= c_3;
      m_2 = m_3;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_2 = 0; j_2 < J_2-1; j_2++){
      for (j_1 = 0; j_1 < J_1; j_1++){
        T[I] += c;
        c *= m_1;
        I += 1;
        T[I] += s;
        s *= m_1;
        I += 1;
      }
      m_2 *= c_2;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_1 = 0; j_1 < J_1-1; j_1++){
      T[I] += c;
      c *= m_1;
      I += 1;
      T[I] += s;
      s *= m_1;
      I += 1;
    }
    
    T[I] += c;
    I += 1;
    T[I] += s;
  }
  
  return T;
}

// [[Rcpp::export()]]

NumericVector ChevTrA(NumericVector T){
  
  int J = T.size(), k, j;
  NumericVector E(J),O(J);
  
  for (j = 0; j < J; j++){
    E(j) = T[j];
    O(j) = T[j];
  }
  
  for (j = 0; j < J; j++){
    T(j) = O[0];
    
    for (k = 0; k < J-j-1; k++){
      E(k) = 2*O[k+1]-E[k];
      O(k) = 2*E[k]-O[k]; 
    }
  }
  
  return T;
}

// [[Rcpp::export()]]

NumericVector ChevTrB(NumericVector T){
  
  int J = T.size(), k, j;
  NumericVector E(J), O(J-1);
  E[0] = T[0];
  
  for (j = 1; j < J; j++){
    E(j) = T[j];
    O(j-1) = T[j];
  }
  
  for (j = 0; j < floor(J/2); j++){
    T(2*j) = E[0];
    T(2*j+1) = O[0];
    
    for (k = 0; k < J-2*j-2; k++){
      E(k) = 2*O[k+1]-E[k];
    }
    
    for (k = 0; k < J-2*j-3; k++){
      O(k) = 2*E[k+1]-O[k];
    }
  }
  
  if (J % 2 == 1){
    T(J-1) = E[0];
  }
  
  return(T);
}

// [[Rcpp::export()]]

NumericVector SineTr(NumericVector T){
  
  int J = T.size(), k, j;
  NumericVector E(J), R(J);
  
  for (k = 1; k < J; k++){
    for (j = k; j < J; j++){
      E(j) = T[j-1]-T[j];
    }
    
    for (j = k; j < J; j++){
      T(j) = E[j];
    }
  }
  
  return(T);
}

// [[Rcpp::export()]]

NumericVector OddNegTr(NumericVector T){
  
  int J = T.size(), j;
  
  for (j = 1; j < J; j = j+2){
    T(j) = -T[j];
  }
  
  return(T);
}

// [[Rcpp::export()]]

NumericVector CosineSumTr(NumericMatrix TT){
  
  int J = TT.nrow(), K = TT.ncol(), j, k;
  NumericVector T(J*K);
  
  for (j = 0; j < J; j++){
    T(j) = TT(j,0);
  }
  
  for (k = 1; k < K; k++){
    for (j = 0; j < J; j++){
      T(k*J+j) = 2*TT(j,k)-T[k*J-j-1];
    }
  }
  
  return(T);
}

// [[Rcpp::export()]]

NumericVector SineSumTr(NumericMatrix TT){
  
  int J = TT.nrow(), K = TT.ncol(), j, k;
  NumericVector T(J*K);
  
  for (j = 0; j < J; j++){
    T(j) = TT(j,0);
  }
  
  for (k = 1; k < K; k++){
    for (j =0; j < J; j++){
      T(k*J+j) = 2*TT(j,k)+T[k*J-j-1];
    }
  }
  
  return(T);
}

// [[Rcpp::export()]]

NumericVector get_quantile_multiplier(double m, int J){
  
  int j;
  NumericVector F(2*J+1);
  double m_j = m, m_inc = 2*m, pi = 3.14159265358979323846;
  F(0) = 0;
  
  for (j = 1; j < 2*J; j += 2){
    F[j] = (2/(pi*j*j))*cos(m_j);
    F[j+1] = (2/(pi*j*j))*sin(m_j);
    m_j += m_inc;
  }
  
  return(F);
}

// [[Rcpp::export()]]

NumericVector get_lowess_multiplier(double x, double h, int J){
  
  int j;
  NumericVector F(2*J+1);
  double x_j = x, x_inc = 2*x, h_j = h, h_inc = 2*h, pi = 3.14159265358979323846;
  F(0) = 0;
  
  for (j = 1; j < 2*J; j += 2){
    F[j] = (4/(pi*j))*cos(x_j)*sin(h_j);
    F[j+1] = (4/(pi*j))*sin(x_j)*sin(h_j);
    x_j += x_inc;
    h_j += h_inc;
  }
  
  return(F);
}

// [[Rcpp::export()]]

NumericMatrix get_loess_powers(NumericVector x, NumericVector y, NumericVector x_0, NumericVector h, int K){
  int n = x.size(), m = x_0.size(), i, j, k;
  double p;
  NumericMatrix P(m,3*K+2);
  
  for (j = 0; j < m; j++){
    for (i = 0; i < n; i++){
      if (x_0[j] - h[j] < x[i] & x_0[j] + h[j] > x[i]){
        p = pow(1-pow(fabs((x[i]-x_0[j])/h[j]),3),3);
        P(j,0) += p;
        for (k = 1; k < 2*K+1; k++){
          p *= x[i]-x_0[j];
          P(j,k) += p;
        }
        p = y[i]*pow(1-pow(fabs((x[i]-x_0[j])/h[j]),3),3);
        P(j,k) += p;
        for (k = 2*K+2; k < 3*K+2; k++){
          p *= x[i]-x_0[j];
          P(j,k) += p;
        }
      }
    }
  }
  return(P);
}

// [[Rcpp::export()]]

NumericVector get_loess_match(NumericVector x, NumericVector x_0, double M){
  int n = x.size(), m = x_0.size(), i, j;
  NumericVector x_p(n);
  
  for (i = 0; i < n; i++){
    for (j = 0; j < m; j++){
      if(x[i] == x_0[j]){
        x_p[i] = x_0[j];
        break;
      } 
    }
    if (j == m){
      x_p[i] = M;
    }
  }
  return(x_p);
}

// [[Rcpp::export()]]

NumericVector get_kdtree_summand_2_parameter(NumericMatrix X, int J_1, int J_2){
  
  int i, j, j_1, j_2, k, m, n, p, J, K, N, P;
  N = X.nrow();
  P = X.ncol();
  J = J_1*J_2;
  K = pow(2*J+1,P);
  NumericVector C((2*J+1)*P), L(pow(2*J+1,P-1)), M(pow(2*J+1,P-1)), T(K);
  double c, s, c_1, s_1, c_2, m_1, m_2;
  
  for (n = 0; n < N; n++){
    i = 0;
    
    for (p = 0; p < P-1; p++){
      C[i] = 1;
      i += 1;
      c_1 = cos(X(n,p));
      c_2 = cos(2*J_1*X(n,p));
      m_1 = c_1*c_1;
      c = c_1;
      s_1 = sqrt(1-m_1);
      s = s_1;
      m_2 = 1;
      
      for (j_2 = 0; j_2 < J_2-1; j_2++){
        for (j_1 = 0; j_1 < J_1; j_1++){
          C[i] = c;
          c *= m_1;
          i += 1;
          C[i] = s;
          s *= m_1;
          i += 1;
        }
        m_2 *= c_2;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      
      for (j_1 = 0; j_1 < J_1-1; j_1++){
        C[i] = c;
        c *= m_1;
        i += 1;
        C[i] = s;
        s *= m_1;
        i += 1;
      }
      C[i] = c;
      i += 1;
      C[i] = s;
      i += 1;
    }
    
    C[i] = 1;
    i += 1;
    c_1 = cos(X(n,p));
    c_2 = cos(2*J_1*X(n,p));
    m_1 = c_1*c_1;
    c = c_1;
    s_1 = sqrt(1-m_1);
    s = s_1;
    m_2 = 1;
    
    for (j_2 = 0; j_2 < J_2-1; j_2++){
      for (j_1 = 0; j_1 < J_1; j_1++){
        C[i] = c;
        c *= m_1;
        i += 1;
        C[i] = s;
        s *= m_1;
        i += 1;
      }
      m_2 *= c_2;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_1 = 0; j_1 < J_1-1; j_1++){
      C[i] = c;
      c *= m_1;
      i += 1;
      C[i] = s;
      s *= m_1;
      i += 1;
    }
    C[i] = c;
    i += 1;
    C[i] = s;
    
    for (m = 0; m < 2*J+1; m++){
      M[m] = C[m];
    }
    K = m;
    
    if (P > 2){
      for (p = 1; p < (P-1); p++){
        i = 0;
        
        for (j = 0; j < 2*J+1; j++){
          for (k = 0; k < K; k++){
            L[i] = M[k]*C[m];
            i += 1;
          }
          m += 1;
        }
        K = i;
        
        for (k = 0; k < K; k++){
          M[k] = L[k];
        }
      }
    }
    
    i = 0;
    
    for (j = 0; j < 2*J+1; j++){
      for (k = 0; k < K; k++){
        T[i] += M[k]*C[m];
        i += 1;
      }
      m += 1;
    }
  }
  return T;
}

// [[Rcpp::export()]]

NumericVector get_kdtree_summand_3_parameter(NumericMatrix X, int J_1, int J_2, int J_3) {
  
  int i, j, j_1, j_2, j_3, k, m, n, p, J, K, N, P;
  N = X.nrow();
  P = X.ncol();
  J = J_1*J_2*J_3;
  K = pow(2*J+1,P);
  NumericVector C((2*J+1)*P), L(pow(2*J+1,P-1)), M(pow(2*J+1,P-1)), T(K);
  double c, s, c_1, s_1, c_2, c_3, m_1, m_2, m_3;
  
  for (n = 0; n < N; n++){
    i = 0;
    
    for (p = 0; p < P-1; p++){
      C[i] = 1;
      i += 1;
      c_1 = cos(X(n,p));
      c_2 = cos(2*J_1*X(n,p));
      c_3 = cos(2*J_2*J_1*X(n,p));
      m_1 = c_1*c_1;
      c = c_1;
      s_1 = sqrt(1-m_1);
      s = s_1;
      m_2 = 1;
      m_3 = 1;
      
      for (j_3 = 0; j_3 < J_3-1; j_3++){
        for (j_2 = 0; j_2 < J_2; j_2++){
          for (j_1 = 0; j_1 < J_1; j_1++){
            C[i] = c;
            c *= m_1;
            i += 1;
            C[i] = s;
            s *= m_1;
            i += 1;
          }
          m_2 *= c_2;
          c = c_1*m_2;
          s = s_1*m_2;
        }
        m_3 *= c_3;
        m_2 = m_3;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      
      for (j_2 = 0; j_2 < J_2-1; j_2++){
        for (j_1 = 0; j_1 < J_1; j_1++){
          C[i] = c;
          c *= m_1;
          i += 1;
          C[i] = s;
          s *= m_1;
          i += 1;
        }
        m_2 *= c_2;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      
      for (j_1 = 0; j_1 < J_1-1; j_1++){
        C[i] = c;
        c *= m_1;
        i += 1;
        C[i] = s;
        s *= m_1;
        i += 1;
      }
      C[i] = c;
      i += 1;
      C[i] = s;
      i += 1;
    }
    
    C[i] = 1;
    i += 1;
    c_1 = cos(X(n,p));
    c_2 = cos(2*J_1*X(n,p));
    c_3 = cos(2*J_2*J_1*X(n,p));
    m_1 = c_1*c_1;
    c = c_1;
    s_1 = sqrt(1-m_1);
    s = s_1;
    m_2 = 1;
    m_3 = 1;
    
    for (j_3 = 0; j_3 < J_3-1; j_3++){
      for (j_2 = 0; j_2 < J_2; j_2++){
        for (j_1 = 0; j_1 < J_1; j_1++){
          C[i] = c;
          c *= m_1;
          i += 1;
          C[i] = s;
          s *= m_1;
          i += 1;
        }
        m_2 *= c_2;
        c = c_1*m_2;
        s = s_1*m_2;
      }
      m_3 *= c_3;
      m_2 = m_3;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_2 = 0; j_2 < J_2-1; j_2++){
      for (j_1 = 0; j_1 < J_1; j_1++){
        C[i] = c;
        c *= m_1;
        i += 1;
        C[i] = s;
        s *= m_1;
        i += 1;
      }
      m_2 *= c_2;
      c = c_1*m_2;
      s = s_1*m_2;
    }
    
    for (j_1 = 0; j_1 < J_1-1; j_1++){
      C[i] = c;
      c *= m_1;
      i += 1;
      C[i] = s;
      s *= m_1;
      i += 1;
    }
    C[i] = c;
    i += 1;
    C[i] = s;
    
    for (m = 0; m < 2*J+1; m++){
      M[m] = C[m];
    }
    K = m;
    
    if (P > 2){
      for (p = 1; p < (P-1); p++){
        i = 0;
        
        for (j = 0; j < 2*J+1; j++){
          for (k = 0; k < K; k++){
            L[i] = M[k]*C[m];
            i += 1;
          }
          m += 1;
        }
        K = i;
        
        for (k = 0; k < K; k++){
          M[k] = L[k];
        }
      }
    }
    
    i = 0;
    
    for (j = 0; j < 2*J+1; j++){
      for (k = 0; k < K; k++){
        T[i] += M[k]*C[m];
        i += 1;
      }
      m += 1;
    }
  }
  return T;
}

// [[Rcpp::export]]

NumericVector get_kdtree_multiplier(float a, float b, int J){
  
  int j;
  NumericVector F(2*J+1);
  double a_j = a, a_inc = 2*a, b_j = b, b_inc = 2*b, pi = 3.141592653589793115998;
  F[0] = 1;
  
  if(a > 0){
    F[0] -= .5;
    for (j = 1; j < 2*J; j += 2){
      F[j] -= 2/(pi*j)*sin(a_j);
      F[j+1] += 2/(pi*j)*cos(a_j);
      a_j += a_inc;
    }
  }
  
  if(b < 1){
    F[0] -= .5;
    for (j = 1; j < 2*J; j += 2){
      F[j] += 2/(pi*j)*sin(b_j);
      F[j+1] -= 2/(pi*j)*cos(b_j);
      b_j += b_inc;
    }
  }
  
  return(F);
}

// [[Rcpp::export]]

NumericVector kdtree_medians2cellcounts(NumericMatrix data, NumericMatrix medians) {
  int n, nn, p, pp, d, dd, i, I;
  n = data.nrow(), p = data.ncol(), d = medians.nrow(), I = medians.ncol();
  NumericVector cell_count(2*I);
  for (nn = 0; nn < n; nn++){
    i = 0;
    for (dd = 0; dd < d; dd++){
      pp = dd % p;
      if (data(nn,pp) < medians(dd,i)){
        i = 2*i;
      } else {
        i = 2*i+1;
      }
    }
    cell_count(i) += 1; 
  }
  return(cell_count);
}