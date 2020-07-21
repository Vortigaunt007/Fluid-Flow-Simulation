#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "Matrix.h"
#include "Vector.h"
#include "Tensor.h"

class Grid
{
    enum typePoint { ANGLE, LEFT, RIGHT, BOTTOM, TOP, CENTER };

    size_t N;
    size_t M;
    double hx, hy;
public:
    Grid();
    Grid(size_t, size_t);

    typePoint type(size_t i, size_t j);
    void setGridStep();

    Matrix initialization_A(Tensor &K);
    Matrix initialization_A_test();
    Vector initialization_b_test();

    Vector initialization_b(Tensor &K, Matrix &A, size_t l_dirichlet_start = 0, size_t l_dirichlet_end = 0, double value_dirichlet = 0, size_t num = 1,
                            size_t l1_start = 0, size_t l1_end = 0, double value1 = 0,
                            size_t l2_start = 0, size_t l2_end = 0, double value2 = 0,
                            size_t l3_start = 0, size_t l3_end = 0, double value3 = 0,
                            size_t l4_start = 0, size_t l4_end = 0, double value4 = 0);



    size_t getPointIndex(size_t, size_t);

    size_t getN() const;
    size_t getM() const;
    double getHx() const;
    double getHy() const;
};

#endif // GRID_H
