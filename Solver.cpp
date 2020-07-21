#include "Solver.h"


Solver::Solver(): eps(0.001)
{

}

Vector Solver::CraigMethod(Grid &G, Matrix &A, Vector &b)
{
    size_t N = G.getN();
    size_t M = G.getM();
    size_t k = N * M;

    Vector x_0(k);
    Vector Ax_k(k);
    Vector x_k(k);
    Vector r_k(k);
    Vector r_km1(k);
    Vector r_1(k);
    Vector q_k(k);
    Vector q_km1(k);
    Vector p_k(k);
    Vector stop(k);
    Vector z_ij(k);

    //b.saveToFile(N, M, "B.txt");

    Matrix A_T(k, k);

    double alpha, beta;


    // Initialization

   // std::cout << x_k << std::endl;

    trans(A, A_T);

    A.fillNonZeroElementsPosition();
    A_T.fillNonZeroElementsPosition();
    x_k = b;
    mult_Matrix_Vector(A, x_k, Ax_k);

   // std::cout << Ax_k << std::endl;
    summAlpha(Ax_k, -1.0, b, r_1);
    r_k = r_1;

    int count_step = 0;

  //  std::cout << r_1 << std::endl;
   // std::cout << r_k << std::endl;

    while ((r_k.normL2() / r_1.normL2()) > eps) {

       // std::cout << r_k.normL2() / r_1.normL2() << std::endl;
std::cout << count_step++ << std::endl;
        //std::cout << A << std::endl;
        //std::cout << x_k << std::endl;

        mult_Matrix_Vector(A, x_k, stop);

      //  std::cout << stop << std::endl;
        summAlpha(stop, -1.0, b, stop);
        if (dotProduct(stop, stop) < 0.000001)
            break;

        alpha = 1.0 / dotProduct(r_k, r_k);
        summAlpha(p_k, alpha, r_k, p_k);
        mult_Matrix_Vector(A_T, p_k, q_k);
        beta = 1.0 / dotProduct(q_k, q_k);
        summAlpha(x_k, -beta, q_k, x_k);
        mult_Matrix_Vector(A, x_k, Ax_k);
        summAlpha(Ax_k, -1.0, b, r_k);
    }


    //std::cout << x_k << std::endl;

    Matrix A_n(k, k);
    Vector z_n(k);
    Vector x_n(k);
    Vector b_n(k);

    A_n.copyFromData(A);
    sqr(A_n);
    x_n = x_k;
    sqr(x_n);
    b_n = b;
    sqr(b_n);

    A_n.fillNonZeroElementsPosition();

    mult_Matrix_Vector(A_n, x_n, z_n);

    double z = dotProduct(z_n, b_n);
    long double eps_craig = 1.0 / 10000000000000000;
    double r_1_craig = 1.0 / dotProduct(r_1, r_1);
    double err = eps_craig * eps_craig * z * r_1_craig;

    while (err <= 1) {
        mult_Matrix_Vector(A, x_k, stop);
        summAlpha(stop, -1.0, b, stop);
        std::cout << count_step++ << std::endl;
        if (dotProduct(stop, stop) < 0.000001)
            break;
        alpha = 1.0 / dotProduct(r_k, r_k);
        summAlpha(p_k, alpha, r_k, p_k);
        mult_Matrix_Vector(A_T, p_k, q_k);
        beta = 1.0 / dotProduct(q_k, q_k);
        summAlpha(x_k, -beta, q_k, x_k);
        mult_Matrix_Vector(A, x_k, Ax_k);
        summAlpha(Ax_k, -1.0, b, r_k);

        r_1_craig += 1.0 / dotProduct(r_k, r_k);
        err = eps_craig * eps_craig * z * r_1_craig;
        std::cout << err << std::endl;
    }

/*
    x_k[0] = 0.0;
    x_k[N-1] = 0.0;
    x_k[(M - 1) * N] = 0.0;
    x_k[M * N - 1] = 0.0;
*/

    x_k[0] = (x_k[1] + x_k[N]) / 2.0;
    x_k[N-1] = (x_k[N - 2] + x_k[2 * N - 1]) / 2.0;
    x_k[(M - 1) * N] = (x_k[(M - 2) * N] + x_k[(M - 1) * N + 1]) / 2.0;
    x_k[M * N - 1] = (x_k[(M - 1) * N - 1] + x_k[M * N - 2]) / 2.0;

    x_k.saveToFile(N, M, "res.txt");
  //  std::cout << x_k << std::endl;

    return x_k;
}
