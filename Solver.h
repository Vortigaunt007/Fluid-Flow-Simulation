#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"

class Solver
{
private:
    double eps;
public:
    Solver();

    double Green_Function(Grid &G, size_t x, size_t y, size_t X, size_t Y);
    double Boundary_Integral(size_t x, size_t y, Grid &G, Vector &b);
    Vector Analytical_solution(Grid &G, Vector &b_0, Vector &b, double index_dir);
    Vector CraigMethod(Grid &G, Matrix &A, Vector &b);
};

#endif // SOLVER_H
