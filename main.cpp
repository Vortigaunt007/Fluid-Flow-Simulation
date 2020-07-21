#include <iostream>
#include "Grid.h"
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"

using namespace std;

int main()
{
    Grid G(50, 50);
    Tensor K(G.getN(), G.getM(), G.getHx(), G.getHy(), 3.0, 1.0);
    //Matrix A = G.initialization_A_test();
    //Vector b = G.initialization_b_test();
    size_t left, right, bottom, top;
    left = G.getN() / 2;
    right = G.getN()/ 2;
    bottom = G.getM() / 2;
    top = G.getM() / 2;
    Matrix A = G.initialization_A(K);
   // Vector b = G.initialization_b(K, A, 0, left - 1, 1.0, 3, left - 2, left + 3, 1.0, right - 2, right + 3, -1.0);
    Vector b = G.initialization_b(K, A, 0, left - 10, 1.0, 3, left - 10, left - 6, 1.0, 0, 0, 0, 0, 0, 0, top - 2, top + 10, -1.0);
   // std::cout << A << std::endl;
   // std::cout << b << std::endl;
    Solver S;
    S.CraigMethod(G, A, b);
    return 0;
}
