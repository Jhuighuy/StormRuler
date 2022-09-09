#pragma once

#include <Storm/Blass/MatrixDense.hpp>
#include <iostream>
#include <vector>


// using namespace std;

/**
 * Класс Nvt-расчета, где лежат все данные и разные шаги алгоритма
 */
class Nvt {
private:

  std::vector<double> z; ///< мольные фракции компонентов двухфазной смеси Ni/N,
                         ///< вычисляются при инициализации класса
  std::vector<double>
      K; ///< константы фазового равновесия («K-values»), находятся итерационно
  std::vector<double>
      x; ///< мольные фракции компонентов в жидкости (они же dzeta_L)
  std::vector<double>
      y; ///< мольные фракции компонентов в газе (они же dzeta_V)
  std::vector<double> mu; ///< вектор фазовых потенциалов
  std::vector<double> n_vec_L;
  std::vector<double> n_vec_V;
  std::vector<std::vector<double>> ksi_vecs;
  std::vector<double> ksi_right;
  std::vector<double> psi_vec;
  double a_L, a_V, b_L,
      b_V; ///< коэффициенты уравнений состояние в жидкости и газе
  double nu; ///< мольная фракция газовой фазы, рассчитывается с помощью z и К
             ///< через задачу Речфорда-Райса
  double T; ///< Температура; постоянная, задается при инициализации
  double N;
  double V;
  double m1, m2, Omega_a,
      Omega_b; ///< Параметры уравнения состояния; задаются при инициализации
  double beta; ///< объёмная доля газовой фазы
  double v_ovrl; ///< V / N
  double n_ovrl; ///< N / V
  int N_phi;
  double sigma;
  double lambda;
  std::vector<double> phi;
  std::vector<double> dWdphi;


  std::vector<double> c_coeff; ///< Вектор коэффициентов c
  int n_iter; ///< Номер итерации. Если 0, то берем н.у. для K формулам, иначе с
              ///< предыдущей итерации

  int n_comp; ///< Число компонент смеси

  const double eps = 1e-10; ///< Точность сравнения с нулем
  const double R = 8.31;    ///< Газовая постоянная


  /**
   * Задача Речфорда-Райса для нахождения nu
   */
  double RR(const std::vector<double>& z, const std::vector<double>& K,
            double nu_init = 0.5);

  /**
   * Функция РР
   */
  double F_RR(const std::vector<double>& z, const std::vector<double>& K,
              double nu);

  /**
   * Производная РР
   */
  double dF_RR(const std::vector<double>& z, const std::vector<double>& K,
               double nu);

  /**
   * Функция нахождения beta
   */
  double Find_Beta(double nu, double a_L, double b_L, double a_V, double b_V,
                   double beta_init = 0.5);

  /**
   * Вычисление давления по уравнению состояния
   */
  double P_EoS(double v, double a, double b);

  /**
   * Производная уравнения состояния по v
   */
  double dP_EoS(double v, double a, double b);

  /**
   * Вычисление химических потенциалов mu_i
   */
  double mu_i(const std::vector<double>& n, int i);

  /**
   * Вычисление f (свободной энергии Гельмгольца)
   */
  double Helm(const std::vector<double>& n);

  /**
   * Вычисление химических потенциалов mu_i excess
   */
  double mu_i_ex(const std::vector<double>& n, int i);

  /**
   * Вычисление мольной фракции компоненты в жидкости
   */
  double x_Kz(double z, double K, double nu);

  /**
   * Вычисление мольной фракции компоненты в газе
   */
  double y_Kz(double z, double K, double nu);

  /**
   * Вычисление правой части для МПИ для нахождения K_i
   */
  void K_eq_Right(const std::vector<double>& z, const std::vector<double>& K,
                  double nu, double beta, std::vector<double>& K_right);

  /**
   * Вычисление распределения фаз: K, nu, beta, потенциалы mu (все кладется в
   * соответствующие поля класса)
   */
  void Find_Phases();

  /**
   * Норма разности 2 векторов
   */
  double norm_vec_diff(const std::vector<double>& v1,
                       const std::vector<double>& v2);

  /**
   * поиск долей компонент и коэффициентов уравнения состояния
   */
  void Find_xyab();

  /**
   * поиск векторов долей жидкой и газовой фаз
   */
  void Find_n_vecs(double& n_L, double& n_V);

  void Find_Psi_Ksi(int N);

  void Fill_D1_InvD2(const std::vector<double>& ksi, std::vector<double>& D1,
                     Storm::DenseMatrix<double>& InvD2);

  void dmu_matrix(const std::vector<double>& n, Storm::DenseMatrix<double>& M);

  void ComputeCahnHillard();

  double Pressure_Eq(const std::vector<double>& n);


public:

  std::vector<std::vector<double>>
      a_coeff; ///< Матрица коэффициентов a уравнения состояния
  std::vector<double> b_coeff; ///< Вектор коэффициентов b уравнения состояния

  Nvt(); ///< конструктор по умолчанию

  void NVT_set_param(std::vector<double>& T_crit, std::vector<double>& P_crit,
                     std::vector<double>& ac_facs,
                     std::vector<std::vector<double>>& k_ij,
                     std::vector<double> N_parts, double V1, double T1,
                     double P_init, double N_mesh);

  double Test_RR(const std::vector<double>& K_test, double nu_init);
  double Test_Beta(double nu_test, double a_L_test, double b_L_test,
                   double a_V_test, double b_V_test, double beta_init);
  void Print_K_values();
  void Test_Find_Phases();
  void Print_CH_params();
  void Test_CH();
};
