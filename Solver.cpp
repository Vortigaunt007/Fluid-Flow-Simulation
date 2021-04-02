#include "Solver.h"


Solver::Solver(): eps(0.001)
{

}

double Solver::Green_Function(Grid &G, size_t i, size_t j, size_t I, size_t J)
{
    double pi = 3.141592653589793;

    double sum = 0.0;
    double coeff; // коэффициент нормализации в квадрате Lx = 1, Ly = 1
    double x = i * G.getHx(), y = j * G.getHy(),
           X = I * G.getHx(), Y = J * G.getHy();

    for (size_t m = 0; m < 50; m++) {
        for (size_t n = 0; n < 50; n++) {
            if (n == 0 && m == 0) continue;
            if (n == 0 || m == 0) coeff = 2.0 / (1.0 * 1.0);
            else coeff = 4.0 / (1.0 * 1.0);

            sum += coeff / (pi * n * pi * n / 1.0 / 1.0 + pi * m * pi * m / 1.0 / 1.0)
                    * cos(pi * n * X / 1.0) * cos(pi * m * Y / 1.0)
                    * cos(pi * n * x / 1.0) * cos(pi * m * y / 1.0);
        }
    }

    return -1.0 * sum;
}

double Solver::Boundary_Integral(size_t x, size_t y, Grid &G, Vector &b)
{
    std::vector <size_t> boundaries_index = G.boundary_points();
    double value = 0.0;
    double G_pm;
    std::vector <size_t> IJ(2);

   // std::cout << boundaries_index.size() << std::endl;

    for (size_t j = 0; j < boundaries_index.size(); j++) {

        size_t index_bound = boundaries_index[j];
        IJ = G.getPointCoordinates(index_bound);
        size_t X = IJ[0], Y = IJ[1];
        G_pm = Green_Function(G, x, y, X, Y);

        value += G_pm * b[index_bound] * G.getHx();

//        std::cout << G_pm << " " <<  b[index_bound] <<std::endl;
    }

    return value;
}

Vector Solver::Solver::Analytical_solution(Grid &G, Vector &b_0, Vector &b, double index_dir)
{
    Vector x_k(G.getN() * G.getM());
    std::vector <size_t> ij(2);
    double value;

    for (size_t i = 0; i < x_k.getSize(); i++){

         ij = G.getPointCoordinates(i);
         size_t x = ij[0], y = ij[1];

         value = Boundary_Integral(x, y, G, b_0);
         x_k[i] = value;
    }

    double const_c;
    const_c = b[index_dir] - x_k[index_dir];

    for (size_t i = 0; i < x_k.getSize(); i++)
        x_k[i] = x_k[i] + const_c;

    return x_k;
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

    trans(A, A_T);

    A.fillNonZeroElementsPosition();
    A_T.fillNonZeroElementsPosition();
    x_k = b;
    mult_Matrix_Vector(A, x_k, Ax_k);

    summAlpha(Ax_k, -1.0, b, r_1);
    r_k = r_1;

    int count_step = 0;

    while ((r_k.normL2() / r_1.normL2()) > eps) {

        if (count_step++ % 1000 == 0)
          std::cout << (r_k.normL2() / r_1.normL2()) << std::endl;

        mult_Matrix_Vector(A, x_k, stop);

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

        if (count_step++ % 1000 == 0)
            std::cout << err << std::endl;
    }

    x_k[0] = (x_k[1] + x_k[N]) / 2.0;
    x_k[N-1] = (x_k[N - 2] + x_k[2 * N - 1]) / 2.0;
    x_k[(M - 1) * N] = (x_k[(M - 2) * N] + x_k[(M - 1) * N + 1]) / 2.0;
    x_k[M * N - 1] = (x_k[(M - 1) * N - 1] + x_k[M * N - 2]) / 2.0;

    x_k.saveToFile(N, M, "res.txt");

    return x_k;
}
