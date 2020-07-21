#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"

class Solver
{
private:
    double eps;
public:
    Solver();
    Vector CraigMethod(Grid &G, Matrix &A, Vector &b);
};

#endif // SOLVER_H
