#include "NVT.h"

double NVT::F_RR(const std::vector<double>& z, const std::vector<double>& K,
                 double nu) {
  double sum = 0.0;
  // for (auto [z_i, K_i])
  for (int i = 0; i < n_comp; ++i) {
    sum += z[i] * (K[i] - 1) / (1 + nu * (K[i] - 1));
  }
  return sum;
};

double NVT::dF_RR(const std::vector<double>& z, const std::vector<double>& K,
                  double nu) {
  double sum = 0.0;
  for (int i = 0; i < n_comp; ++i) {
    sum += z[i] * (K[i] - 1) * (K[i] - 1) /
           ((1 + nu * (K[i] - 1)) * (1 + nu * (K[i] - 1)));
  }
  return -sum;
};

double NVT::RR(const std::vector<double>& z, const std::vector<double>& K,
               double nu_init) {
  //Сначала проверка наличия корня, если его нет, то выдаем 0
  double const F_0 = F_RR(z, K, 0);
  double const F_1 = F_RR(z, K, 1);
  if (F_0 * F_1 > 0.0) {
    if (std::abs(F_0) < eps) return 0.0;
    if (std::abs(F_1) < eps) return 1.0;
    return 0.0;
  }

  double nu = nu_init;
  double nu_tmp;
  double F = F_RR(z, K, nu);
  while (std::abs(F) > eps) {
    nu_tmp = nu - F / dF_RR(z, K, nu);
    //Если nu вышло за границы, то надо вместо Ньютона сделать шаг по бисекции
    if ((nu_tmp < 1) && (nu_tmp > 0)) {
      nu = nu_tmp;
    } else {
      if (nu_tmp > 1) {
        nu = 0.5 * (nu + 1);
      } else {
        nu = 0.5 * nu;
      }
    }
    F = F_RR(z, K, nu);
  }

  return nu;
};

double NVT::Find_Beta(double nu, double a_L, double b_L, double a_V, double b_V,
                      double beta_init) {
  //Простые случаи
  if (nu <= eps) return 0.0;
  if (1 - nu <= eps) return 1.0;

  int step_max = 200000;
  int step = 0;

  double beta = beta_init;
  double beta_tmp;
  double F = P_EoS(beta * v_ovrl / nu, a_V, b_V) -
             P_EoS((1 - beta) * v_ovrl / (1 - nu), a_L, b_L);
  while (std::abs(F) > eps && (step < step_max)) {
    /*beta_tmp = beta - F / (v_ovrl / nu *  dP_EoS(beta * v_ovrl / nu, a_V, b_V)
       + v_ovrl / (1 - nu) * dP_EoS((1 - beta) * v_ovrl / (1 - nu), a_L,
       b_L));*/
    beta_tmp =
        beta -
        F / (1.0 / nu * dP_EoS(beta * v_ovrl / nu, a_V, b_V) +
             1.0 / (1 - nu) * dP_EoS((1 - beta) * v_ovrl / (1 - nu), a_L, b_L));
    //Если beta вышло за границы beta = 1 и beta = 0, то надо вместо Ньютона
    //сделать шаг по бисекции
    if (beta_tmp > 1) {
      beta = 0.5 * (beta + 1);
    } else {
      if (beta_tmp < 0) {
        beta = 0.5 * beta;
      } else {
        beta = beta_tmp;
      }
    }
    F = P_EoS(beta * v_ovrl / nu, a_V, b_V) -
        P_EoS((1 - beta) * v_ovrl / (1 - nu), a_L, b_L);

    step++;
  }

  return beta;
};

double NVT::P_EoS(double v, double a, double b) {
  return R * T / (v - b) - a / (v + m1 * b) / (v + m2 * b);
};

double NVT::dP_EoS(double v, double a, double b) {
  return -R * T / (v - b) / (v - b) + a / (b * (m2 - m1)) *
                                          (1.0 / (v + m1 * b) / (v + m1 * b) -
                                           1.0 / (v + m2 * b) / (v + m2 * b));
};

// double NVT::mu_i(const vector<double>& n, int i, double n_sum, double a,
// double b) {
//     double sum = 0.0;
//     for (int j = 0; j < n_comp; ++j) {
//         sum += a_coeff[i][j] * n[j];
//     }
//     double res = (2 * sum / n_sum - a * b_coeff[i] / b) / b / (m2 - m1) +
//         R * T * (log(n[i] / (1 - b * n_sum)) + b_coeff[i] * n_sum / (1 - b *
//         n_sum))
//         - a * n_sum * b_coeff[i] / b / (1 + m1 * b * n_sum) / (1 + m2 * b *
//         n_sum);
//     return (2 * sum / n_sum - a * b_coeff[i] / b) / b / (m2 - m1) +
//         R * T * (log(n[i] / (1 - b * n_sum)) + b_coeff[i] * n_sum / (1 - b *
//         n_sum))
//         - a * n_sum * b_coeff[i] / b / (1 + m1 * b * n_sum) / (1 + m2 * b *
//         n_sum);
// }

double NVT::mu_i(const std::vector<double>& n, int i) {
  double dn = 1e-3;

  std::vector<double> n_dn_right;
  std::vector<double> n_dn_left;
  n_dn_right = n;
  n_dn_left = n;
  double n_sum_right, n_sum_left;
  n_dn_right[i] = n[i] * (1.0 + dn);
  n_dn_left[i] = n[i] * (1.0 - dn);
  // n_sum_right = n_sum + n[i] * dn;
  // n_sum_left = n_sum - n[i] * dn;
  // double tmp = (Helm(n_dn_right, n_sum_right, a, b) - Helm(n_dn_left,
  // n_sum_left, a, b)) / (n_dn_right[i] - n_dn_left[i]);
  double tmp =
      (Helm(n_dn_right) - Helm(n_dn_left)) / (n_dn_right[i] - n_dn_left[i]);
  // return (Helm(n_dn_right, n_sum_right, a, b) - Helm(n_dn_left, n_sum_left,
  // a, b)) / (n_dn_right[i] - n_dn_left[i]);
  return tmp;
}

void NVT::dmu_matrix(const std::vector<double>& n,
                     Storm::DenseMatrix<double>& M) {
  double dn = 1e-3;
  for (int i = 0; i < n_comp; ++i) {
    std::vector<double> n_dn_right(n);
    std::vector<double> n_dn_left(n);
    n_dn_right[i] = n[i] * (1.0 + dn);
    n_dn_left[i] = n[i] * (1.0 - dn);
    double tmp;
    for (int j = 0; j < n_comp; ++j) {
      tmp = (mu_i(n_dn_right, j) - mu_i(n_dn_left, j)) /
            (n_dn_right[i] - n_dn_left[i]) / c_coeff[i] / c_coeff[j];
      M(i, j) = (mu_i(n_dn_right, j) - mu_i(n_dn_left, j)) /
                (n_dn_right[i] - n_dn_left[i]) / c_coeff[i] / c_coeff[j];
    }
  }
}

double NVT::mu_i_ex(const std::vector<double>& n, int i) {
  double tmp = mu_i(n, i);
  // return mu_i(n, i, n_sum, a, b) - R * T * log(n[i]);
  return tmp - R * T * std::log(n[i]);
}

double NVT::x_Kz(double z, double K, double nu) {
  return z / (1.0 + nu * (K - 1.0));
}

double NVT::y_Kz(double z, double K, double nu) {
  return z * K / (1.0 + nu * (K - 1.0));
}

void NVT::K_eq_Right(const std::vector<double>& z, const std::vector<double>& K,
                     double nu, double beta, std::vector<double>& K_right) {
  double n_L = 0.0;
  double n_V = 0.0;

  //Вообще обработку чистого газа и чистой жидкости надо проводить отдельно и не
  //залезать в эту функцию if ((nu <= eps) || (beta <= eps)) {
  //    n_V = 0.0;
  //    n_L = n_ovrl;
  //    //нет газа => все K_i = 0
  //    return;
  //    }
  // else {
  //    if ((1.0 - nu <= eps) || (1.0 - beta <= eps)) {
  //        n_V = n_ovrl;
  //        n_L = 0.0;
  //        //Вообще в таком случае все K_i = inf
  //    }
  //    else {
  //        n_L = n_ovrl * (1.0 - nu) / (1.0 - beta);
  //        n_V = n_ovrl * nu / beta;
  //    }
  //}

  ////vector<double> n_vec_L, n_vec_V;
  // for (int i = 0; i < n_comp; ++i) {
  //     n_vec_L[i] = n_L * x[i];
  //     n_vec_V[i] = n_V * y[i];
  // }

  Find_n_vecs(n_L, n_V);

  //Строим правую часть
  double ln_n = std::log(n_L / n_V);
  double mu_ex_L, mu_ex_V;
  for (int i = 0; i < n_comp; ++i) {
    // mu_ex_V = mu_i_ex(n_vec_V, i, n_V, a_V, b_V);
    // mu_ex_L = mu_i_ex(n_vec_L, i, n_L, a_L, b_L);
    mu_ex_V = mu_i_ex(n_vec_V, i);
    mu_ex_L = mu_i_ex(n_vec_L, i);
    // K_right[i] = ln_n + (mu_i_ex(n_vec_L, i, n_L, a_L, b_L) -
    // mu_i_ex(n_vec_V, i, n_V, a_V, b_V)) / (R * T);
    K_right[i] = ln_n + (mu_ex_L - mu_ex_V) / (R * T);
  }

  // n_vec_L.clear();
  // n_vec_V.clear();

  return;
}

void NVT::Find_Phases() {
  // if (n_iter == 0) {
  //     //Находим н.п. для K[i] по формулам и отдельно для них nu и beta
  // }
  // else {
  //     //Н.п. уже лежат в нужных переменных

  //}

  //Первичный поиск nu, nuta, a, b, x, y
  /*nu = RR(z, K, nu);
  Find_xyab();
  beta = Find_Beta(nu, a_L, b_L, a_V, b_V, beta);*/

  std::vector<double> K_next(K);
  std::vector<double> K_right(n_comp, 0.0);
  // K_eq_Right(z, K, nu, beta, K_right);
  // for (size_t i{0}; i < n_comp; ++i) {
  //     K_next[i] = exp(K_right[i]);
  // }
  int max_iter = 50;
  int iter = 0;
  double norm = 10 * eps;
  double n_L, n_V;
  //Критерий останова - норма разности приближений <= eps
  while (norm > eps && iter < max_iter) {
    std::swap(K, K_next);
    //Находим nu для новых K_i
    nu = RR(z, K, nu);
    // nu = 0.750845771005761;
    //Находим x, y, a, b
    Find_xyab();
    //Находим beta
    beta = Find_Beta(nu, a_L, b_L, a_V, b_V, beta);
    // beta = 0.789653053634341;

    //Находим правую часть уравнения для K_i
    K_eq_Right(z, K, nu, beta, K_right);
    //Новые K_i
    for (int i = 0; i < n_comp; ++i) {
      K_next[i] = std::exp(K_right[i]);
    }
    iter++;
    norm = norm_vec_diff(K_next, K);
  }
  // K найдены, осталось найти для этого K nu и beta (для следующего шага)
  std::swap(K_next, K);
  nu = RR(z, K, nu);
  Find_xyab();
  beta = Find_Beta(nu, a_L, b_L, a_V, b_V, beta);
  Find_n_vecs(n_L, n_V);

  //Считаем итоговые фазовые потенциалы
  for (int i = 0; i < n_comp; ++i) {
    mu[i] = 0.5 * (mu_i(n_vec_L, i) + mu_i(n_vec_V, i));
  }
}

double NVT::norm_vec_diff(const std::vector<double>& v1,
                          const std::vector<double>& v2) {
  double norm = 0.0;
  int size = v1.size();
  for (int i = 0; i < size; ++i) {
    norm += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }
  return std::sqrt(norm);
}

void NVT::Find_xyab() {
  //По значениям K_i, z_i ищутся x_i, y_i, a_L, a_V, b_L, b_V
  for (int i = 0; i < n_comp; ++i) {
    x[i] = x_Kz(z[i], K[i], nu);
    y[i] = y_Kz(z[i], K[i], nu);
  }
  a_L = 0.0;
  a_V = 0.0;
  b_L = 0.0;
  b_V = 0.0;
  for (int i = 0; i < n_comp; ++i) {
    b_L += x[i] * b_coeff[i];
    b_V += y[i] * b_coeff[i];
    for (int j = 0; j < n_comp; ++j) {
      a_L += x[i] * x[j] * a_coeff[i][j];
      a_V += y[i] * y[j] * a_coeff[i][j];
    }
  }
  return;
}

double NVT::Helm(const std::vector<double>& n) {
  double sum = 0.0, n_sum = 0.0;
  double a = 0, b = 0;
  // vector<double> n_new{ 7.549988305190146e3, 0.148737034710701e3,
  // 0.001747643370334e3, 0.000038110015138e3, 0.000000519247902e3 }; a =
  // 0.214008571376308; b = 2.769115505782897e-5;
  for (int i = 0; i < n_comp; ++i) {
    n_sum += n[i];
    sum += n[i] * (std::log(n[i]) - 1.0);
  }
  for (int i = 0; i < n_comp; ++i) {
    b += b_coeff[i] * n[i];
    for (int j = 0; j < n_comp; ++j) {
      a += a_coeff[i][j] * n[i] * n[j];
    }
  }
  b /= n_sum;
  a /= (n_sum * n_sum);
  double tmp = R * T * sum - n_sum * R * T * std::log(1.0 - b * n_sum) +
               std::log((1.0 + m1 * b * n_sum) / (1.0 + m2 * b * n_sum)) * a *
                   n_sum / b / (m2 - m1);
  return R * T * sum - n_sum * R * T * std::log(1.0 - b * n_sum) +
         std::log((1.0 + m1 * b * n_sum) / (1.0 + m2 * b * n_sum)) * a * n_sum /
             b / (m2 - m1);
}

// double NVT::Calc_P_E(const vector<double>& n, double n_sum, double a, double
// b) {
//     double sum = 0.0;
//     //Нужные mu уже лежат в соответствующем поле
//     for (int i = 0; i < n_comp; ++i) {
//         sum += n[i] * mu[i];
//     }
//     return -f(n, n_sum, a, b) + sum;
// }
//
// void NVT::Build_Ksi_Mesh() {
//
// }

void NVT::Find_n_vecs(double& n_L, double& n_V) {
  // double n_V, n_L;
  if ((nu <= eps) || (beta <= eps)) {
    n_V = 0.0;
    n_L = n_ovrl;
    //нет газа => все K_i = 0
    return;
  } else {
    if ((1.0 - nu <= eps) || (1.0 - beta <= eps)) {
      n_V = n_ovrl;
      n_L = 0.0;
      //Вообще в таком случае все K_i = inf
    } else {
      n_L = n_ovrl * (1.0 - nu) / (1.0 - beta);
      n_V = n_ovrl * nu / beta;
    }
  }

  // vector<double> n_vec_L, n_vec_V;
  for (int i = 0; i < n_comp; ++i) {
    n_vec_L[i] = n_L * x[i];
    n_vec_V[i] = n_V * y[i];
  }
}

NVT::NVT_set_param(std::vector<double>& T_crit, std::vector<double>& P_crit,
                   std::vector<double>& ac_facs,
                   std::vector<std::vector<double>>& k_ij,
                   std::vector<double> N_parts, double V1, double T1,
                   double P_init, double N_mesh) {
  T = T1;
  V = V1;
  N = 0.0;
  N_phi = N_mesh;
  n_comp = N_parts.size();
  for (int i = 0; i < n_comp; ++i) {
    N += N_parts[i];
  }

  m1 = 1.0 - sqrt(2.0);
  m2 = 1.0 + sqrt(2.0);
  Omega_a = 0.45724;
  Omega_b = 0.0778;
  std::vector<double> a_ind(n_comp);
  std::vector<double> kappa(n_comp);
  // vector<double> c_coeff(n_comp);
  b_coeff = std::vector<double>(n_comp);
  a_coeff = std::vector<std::vector<double>>();
  c_coeff = std::vector<double>(n_comp);
  double A, B;
  for (int i = 0; i < n_comp; ++i) {
    if (ac_facs[i] < 0.49) {
      kappa[i] =
          0.37464 + 1.54226 * ac_facs[i] - 0.26992 * ac_facs[i] * ac_facs[i];
    } else {
      kappa[i] = 0.37964 + 1.48503 * ac_facs[i] -
                 0.16442 * ac_facs[i] * ac_facs[i] +
                 0.01667 * ac_facs[i] * ac_facs[i] * ac_facs[i];
    }
    // kappa[i] =
    //     0.37464 + 1.54226 * ac_facs[i] - 0.26992 * ac_facs[i] * ac_facs[i];
    a_ind[i] = (1.0 + kappa[i] * (1.0 - std::sqrt(T / T_crit[i])));
    a_ind[i] = a_ind[i] * a_ind[i] * Omega_a * R * R * T_crit[i] * T_crit[i] /
               P_crit[i];
    b_coeff[i] = Omega_b * R * T_crit[i] / P_crit[i];
    A = -1e-16 / (1.2326 + 1.3757 * ac_facs[i]);
    B = 1e-16 / (0.9051 + 1.5410 * ac_facs[i]);
    c_coeff[i] = std::sqrt(a_ind[i] * std::pow(b_coeff[i], 2.0 / 3.0) *
                           (A * (1.0 - T / T_crit[i]) + B));
    // c_coeff[i] = a_ind[i] * pow(b_coeff[i], 2.0 / 3.0) * (A * (1.0 - T /
    // T_crit[i]) + B);
  }
  for (int i = 0; i < n_comp; ++i) {
    std::vector<double> tmp(n_comp);
    for (int j = 0; j < n_comp; ++j) {
      tmp[j] = std::sqrt(a_ind[i] * a_ind[j]) * (1.0 - k_ij[i][j]);
    }
    a_coeff.push_back(tmp);
  }

  n_iter = 0;

  v_ovrl = V / N;
  n_ovrl = 1.0 / v_ovrl;
  nu = 0.75;
  beta = 0.75;

  K = std::vector<double>(n_comp);
  z = std::vector<double>(n_comp);
  x = std::vector<double>(n_comp);
  y = std::vector<double>(n_comp);
  mu = std::vector<double>(n_comp);
  n_vec_L = std::vector<double>(n_comp);
  n_vec_V = std::vector<double>(n_comp);
  ksi_right = std::vector<double>(N_phi);
  psi_vec = std::vector<double>(N_phi);
  ksi_vecs = std::vector<std::vector<double>>();
  for (int i = 0; i < N_phi; ++i) {
    std::vector<double> tmp(n_comp);
    ksi_vecs.push_back(tmp);
  }
  phi = std::vector<double>(N_phi);
  for (int i = 0; i < N_phi; ++i) {
    phi[i] = i * (1.0 / (N_phi - 1));
  }
  dWdphi = std::vector<double>(N_phi);

  //н у для K_i
  for (int i = 0; i < n_comp; ++i) {
    // K[i] = exp(5.37 * (1.0 + ac_facs[i]) * (1.0 - T / T_crit[i]) - log(P_init
    // / P_crit[i]));
    K[i] = exp(5.37 * (1.0 + ac_facs[i]) * (1.0 - T_crit[i] / T) -
               log(P_init / P_crit[i]));
    z[i] = N_parts[i] / N;
  }
  // nu = RR(z, K, nu);
  Find_xyab();
  // beta = Find_Beta(nu, a_L, b_L, a_V, b_V, beta);
}


NVT::NVT() {}

double NVT::Test_RR(const std::vector<double>& K_test, double nu_init) {
  return RR(z, K_test, nu_init);
}

double NVT::Test_Beta(double nu_test, double a_L_test, double b_L_test,
                      double a_V_test, double b_V_test, double beta_init) {
  return Find_Beta(nu_test, a_L_test, b_L_test, a_V_test, b_V_test, beta_init);
}

void NVT::Print_K_values() {
  std::cout << "K-values: ";
  for (auto el : K) {
    std::cout << el << ' ';
  }
  std::cout << endl;
}

void NVT::Test_Find_Phases() {
  Find_Phases();
}

void NVT::Find_Psi_Ksi(int N) {
  // phi in [0,1]; phi[i] = i / (N - 1)

  std::vector<double> eta(n_comp);
  double tmp;
  double ksi_V = 0.0, ksi_L = 0.0;
  for (int i = 0; i < n_comp; ++i) {
    ksi_V += c_coeff[i] * n_vec_V[i];
    ksi_L += c_coeff[i] * n_vec_L[i];
    eta[i] = mu[i] / c_coeff[i];
  }
  tmp = (ksi_L - ksi_V) / (N - 1);
  for (int i = 0; i < N; ++i) {
    ksi_right[i] = ksi_V + tmp * i;
  }

  //Память на ksi_vecs выделена и в случае чего может быть изменена
  for (int i = 0; i < n_comp; ++i) {
    ksi_vecs[0][i] = c_coeff[i] * n_vec_V[i];
    ksi_vecs[N - 1][i] = c_coeff[i] * n_vec_L[i];
  }
  psi_vec[0] = 0.0;
  psi_vec[N - 1] = 0.0;

  Storm::DenseMatrix<double> Theta(n_comp, n_comp);
  std::vector<double> df(n_comp);
  std::vector<double> theta_i(n_comp);
  double theta_sum;
  std::vector<double> ksi_new(n_comp);
  double psi_new;

  int max_iter = 10;

  //Цикл по N - 2 значениям ksi_righr
  for (int i = 1; i < N - 1; ++i) {
    Fill_D1_InvD2(ksi_vecs[i - 1], df, Theta);
    theta_sum = 0.0;
    // fill(theta_i.begin(), theta_i.end(), 0.0);
    for (int j = 0; j < n_comp; ++j) {
      theta_i[j] = 0.0;
      for (int k = 0; k < n_comp; ++k) {
        theta_i[j] += Theta(j, k);
      }
      theta_sum += theta_i[j];
    }

    // std::cout << Theta << std::endl;

    for (int j = 0; j < n_comp; ++j) {
      ksi_vecs[i][j] =
          ksi_vecs[i - 1][j] +
          theta_i[j] / theta_sum * (ksi_right[i] - ksi_right[i - 1]);
    }

    //Запускаем итерационный процесс
    int iter = 0;
    double norm = 10 * eps;
    while ((norm > eps) && (iter < max_iter)) {
      Fill_D1_InvD2(ksi_vecs[i], df, Theta);
      theta_sum = 0.0;
      for (int j = 0; j < n_comp; ++j) {
        theta_i[j] = 0.0;
        for (int k = 0; k < n_comp; ++k) {
          theta_i[j] += Theta(j, k);
        }
        theta_sum += theta_i[j];
      }
      psi_new = 0.0;
      for (int j = 0; j < n_comp; ++j) {
        ksi_new[j] = ksi_vecs[i][j];
        for (int k = 0; k < n_comp; ++k) {
          ksi_new[j] += (theta_i[j] * theta_i[k] / theta_sum - Theta(j, k)) *
                        (df[k] - eta[k]);
        }
        psi_new += theta_i[j] * (df[j] - eta[j]);
      }
      psi_new /= theta_sum;

      norm = norm_vec_diff(ksi_new, ksi_vecs[i]);
      iter++;
      std::swap(ksi_new, ksi_vecs[i]);
    }
    // ksi уже лежат в ksi_vecs
    ;
    psi_vec[i] = psi_new;
  }
}

void NVT::Fill_D1_InvD2(const std::vector<double>& ksi, std::vector<double>& D1,
                        Storm::DenseMatrix<double>& InvD2) {
  std::vector<double> n_loc(ksi);
  for (int i = 0; i < n_comp; ++i) {
    n_loc[i] /= c_coeff[i];
  }
  for (int i = 0; i < n_comp; ++i) {
    D1[i] = mu_i(n_loc, i) / c_coeff[i];
  }
  Storm::DenseMatrix<double> D2(n_comp, n_comp);
  dmu_matrix(n_loc, D2);
  Storm::inplace_inverse_lu(InvD2, D2);

  return;
}

void NVT::ComputeCahnHillard() {
  //Посчитаем равновесное давление по хим потенциалам и фазам
  double Pr = 0.5 * (Pressure_Eq(n_vec_L) + Pressure_Eq(n_vec_V));

  //Считаем интеграл sigma по формуле трапеций
  std::vector<double> n_loc(n_comp);
  double I;
  double I_pred = 0.0;
  sigma = 0.0;
  for (int i = 1; i < N_phi; ++i) {
    n_loc = ksi_vecs[i];
    for (int j = 0; j < n_comp; ++j) {
      n_loc[j] /= c_coeff[j];
    }
    I = Helm(n_loc) + Pr;
    for (int j = 0; j < n_comp; ++j) {
      I -= n_loc[j] * mu[j];
    }
    if (I < 0) I = 0;
    else
      I = sqrt(I);

    sigma += 0.5 * (I + I_pred) * (ksi_right[i] - ksi_right[i - 1]);

    I_pred = I;
  }
  sigma *= sqrt(2.0);

  //считаем lambda
  lambda = (ksi_right[N_phi - 1] - ksi_right[0]) *
           (ksi_right[N_phi - 1] - ksi_right[0]) / sigma;

  //Считаем dW/dphi
  double tmp = lambda * std::sqrt(lambda / sigma);
  for (int i = 0; i < N_phi; ++i) {
    dWdphi[i] = tmp * psi_vec[i];
  }
}

double NVT::Pressure_Eq(const vector<double>& n) {
  double res = -Helm(n);
  for (int i = 0; i < n_comp; ++i) {
    res += n[i] * mu[i];
  }
  return res;
}

void NVT::Print_CH_params() {
  std::cout << "Sigma = " << sigma << std::endl;
  std::cout << "Lambda = " << lambda << std::endl;
  std::cout << "dWdPhi: " << std::endl;
  for (auto el : dWdphi) {
    std::cout << el << std::endl;
  }
}

void NVT::Test_CH() {
  Find_Psi_Ksi(N_phi);
  ComputeCahnHillard();
}
