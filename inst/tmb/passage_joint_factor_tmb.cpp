#include <TMB.hpp>

template<class Type>
Type passage_corr_tmb(Type dist, Type phi, int kernel_code) {
  Type x = dist / phi;
  if (kernel_code == 1) {
    return exp(-x);
  }
  if (kernel_code == 2) {
    Type z = sqrt(Type(3.0)) * x;
    return (Type(1.0) + z) * exp(-z);
  }
  if (kernel_code == 3) {
    Type z = sqrt(Type(5.0)) * x;
    return (Type(1.0) + z + z * z / Type(3.0)) * exp(-z);
  }
  return exp(-Type(0.5) * x * x);
}

template<class Type>
Type euclidean_dist_tmb(matrix<Type> coords, int i, int j) {
  Type out = Type(0.0);
  int d = coords.cols();
  for (int h = 0; h < d; h++) {
    Type diff = coords(i, h) - coords(j, h);
    out += diff * diff;
  }
  return sqrt(out);
}

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_MATRIX(Y);
  DATA_MATRIX(coords);
  DATA_MATRIX(X);
  DATA_IMATRIX(NN);
  DATA_SCALAR(min_tau2);
  DATA_INTEGER(kernel_code);

  PARAMETER_MATRIX(B);
  PARAMETER_VECTOR(log_A_diag);
  PARAMETER_VECTOR(A_free);
  PARAMETER_VECTOR(log_phi);
  PARAMETER_VECTOR(log_tau2);
  PARAMETER_MATRIX(U);

  int n = Y.rows();
  int p = Y.cols();
  int q = X.cols();
  int K = log_phi.size();
  int m = NN.cols();

  matrix<Type> A(p, K);
  A.setZero();
  int idx = 0;
  for (int k = 0; k < K; k++) {
    for (int j = k; j < p; j++) {
      if (j == k) {
        A(j, k) = exp(log_A_diag(k));
      } else {
        A(j, k) = A_free(idx);
        idx++;
      }
    }
  }

  vector<Type> phi(K);
  for (int k = 0; k < K; k++) phi(k) = exp(log_phi(k));

  vector<Type> tau2(p);
  for (int j = 0; j < p; j++) tau2(j) = exp(log_tau2(j)) + min_tau2;

  Type nll = Type(0.0);
  Type jitter = Type(1e-8);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      Type mu = Type(0.0);
      for (int h = 0; h < q; h++) mu += X(i, h) * B(h, j);
      for (int k = 0; k < K; k++) mu += A(j, k) * U(i, k);
      nll -= dnorm(Y(i, j), mu, sqrt(tau2(j)), true);
    }
  }

  for (int k = 0; k < K; k++) {
    for (int i = 0; i < n; i++) {
      int r = 0;
      for (int a = 0; a < m; a++) {
        if (NN(i, a) >= 0) r++;
      }

      if (r == 0) {
        nll -= dnorm(U(i, k), Type(0.0), Type(1.0), true);
      } else {
        std::vector<int> nb(r);
        int pos = 0;
        for (int a = 0; a < m; a++) {
          if (NN(i, a) >= 0) {
            nb[pos] = NN(i, a);
            pos++;
          }
        }

        matrix<Type> C(r, r);
        vector<Type> c(r);
        vector<Type> u_nb(r);
        matrix<Type> rhs_c(r, 1);
        matrix<Type> rhs_u(r, 1);

        for (int a = 0; a < r; a++) {
          int ia = nb[a];
          Type dia = euclidean_dist_tmb(coords, i, ia);
          c(a) = passage_corr_tmb(dia, phi(k), kernel_code);
          u_nb(a) = U(ia, k);
          rhs_c(a, 0) = c(a);
          rhs_u(a, 0) = u_nb(a);

          for (int b = 0; b < r; b++) {
            int ib = nb[b];
            Type dab = euclidean_dist_tmb(coords, ia, ib);
            C(a, b) = passage_corr_tmb(dab, phi(k), kernel_code);
          }
          C(a, a) += jitter;
        }

        matrix<Type> Cinv_c = C.ldlt().solve(rhs_c);
        matrix<Type> Cinv_u = C.ldlt().solve(rhs_u);

        Type cond_mean = Type(0.0);
        Type cond_var_reduction = Type(0.0);
        for (int a = 0; a < r; a++) {
          cond_mean += c(a) * Cinv_u(a, 0);
          cond_var_reduction += c(a) * Cinv_c(a, 0);
        }

        Type cond_var = Type(1.0) - cond_var_reduction;
        if (cond_var < jitter) cond_var = jitter;
        nll -= dnorm(U(i, k), cond_mean, sqrt(cond_var), true);
      }
    }
  }

  REPORT(A);
  REPORT(phi);
  REPORT(tau2);
  REPORT(B);
  REPORT(U);

  ADREPORT(A);
  ADREPORT(phi);
  ADREPORT(tau2);

  return nll;
}
