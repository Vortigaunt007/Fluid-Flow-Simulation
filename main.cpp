#include <iostream>
#include "Grid.h"
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include <iomanip>

using namespace std;

int main()
{
    Grid G(60, 60);
    Tensor K(G.getN(), G.getM(), G.getHx(), G.getHy(), 1.0, 1.0);

    size_t left, right, bottom, top;
    left = G.getN() / 2;
    right = G.getN()/ 2;
    bottom = G.getM() / 2;
    top = G.getM() / 2;
    size_t index_dirichlet = G.getPointIndex(bottom - 1, 0);
    Matrix A = G.initialization_A(K);
    Vector b = G.initialization_b(K, A, index_dirichlet, 1.0, 0, 0, 0, 0, 0, 0, bottom - 5, bottom + 5, 0.0, top - 5, top + 5, 0.0);
    Vector b_0 = G.initialization_b(K, A, 0, 0.0, 0, 0, 0, 0, 0, 0, bottom - 5, bottom + 5, 0.0, top - 5, top + 5, 0.0);

    Solver S;
    Vector x_k(G.getM() * G.getN()), x(G.getM() * G.getN());
    x_k = S.CraigMethod(G, A, b);
    x_k.saveToFile(G.getN(), G.getM(), "method_x.txt");

    // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    x = S.Analytical_solution(G, b_0, b, index_dirichlet);
    x.saveToFile(G.getN(), G.getM(), "analyt_x.txt");

    Vector x_wo_bound = G.Vector_wo_boundary_nodes(x);
    Vector x_k_wo_bound = G.Vector_wo_boundary_nodes(x_k);
    x_wo_bound.saveToFile(G.getN() - 2, G.getM() - 2, "analyt_x2.txt");
    x_k_wo_bound.saveToFile(G.getN() - 2, G.getM() - 2, "method_x2.txt");


    summAlpha(x_wo_bound, -1.0, x_k_wo_bound, x_wo_bound);
    std::cout << "Delta    " << x_wo_bound.normL2() * G.getHx() << std::endl;
    summAlpha(x, -1.0, x_k, x);
    std::cout << "Delta    " << x.normL2() * G.getHx() << std::endl; // 0.00349
/*
    std::vector <double> x_delta;

    summAlpha(x, -1.0, x_k, x);
    x.saveToFile(G.getN(), G.getM(), "delta_x.txt");


    std::vector <size_t> index_bound = G.boundary_points();
    double val_delta;
    for (size_t i = 0; i < x.getSize(); i++) {
        bool flag = 0;
        for (size_t j = 0; j < index_bound.size(); j++) // проверяем, является ли точка граничной
            if (i == index_bound[j]) flag = 1;
        if(!flag) {
            double res = abs(x[i] - x_k[i]);
            val_delta += sqrt(res * res * G.getHx() * G.getHy());
        }
    }

    std::cout << "Delta  =  " << val_delta << std::endl;
/*
    Vector x_D(x_delta.size());
    for (size_t i = 0; i < x_D.getSize(); i++)
        x_D[i] = x_delta[i];

   // x_D.saveToFile(G.getN() - 2, G.getM() -2, "delta_x.txt");

    std::cout << "Delta    " << x_D.normL2() * G.getHx() << std::endl;
*/
    return 0;
}
