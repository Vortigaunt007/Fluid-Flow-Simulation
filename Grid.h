#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "Matrix.h"
#include "Vector.h"
#include "Tensor.h"

class Grid
{
    enum typePoint { ANGLE, LEFT, RIGHT, BOTTOM, TOP, CENTER };

    size_t N; // Oy
    size_t M; // Ox
    double hx, hy;
public:
    Grid();
    Grid(size_t, size_t);

    typePoint type(size_t i, size_t j);
    void setGridStep();

    double coef_a_x(size_t i, size_t j, Tensor &K);
    double coef_a_y(size_t i, size_t j, Tensor &K);

    Matrix initialization_A(Tensor &K);
    Matrix initialization_A_2(Tensor &K);
    Matrix initialization_A_test();
    Vector initialization_b_test();

    Vector initialization_b(Tensor &K, Matrix &A, size_t index_dirichlet = 0, double value_dirichlet = 0,
                            size_t l1_start = 0, size_t l1_end = 0, double value1 = 0,
                            size_t l2_start = 0, size_t l2_end = 0, double value2 = 0,
                            size_t l3_start = 0, size_t l3_end = 0, double value3 = 0,
                            size_t l4_start = 0, size_t l4_end = 0, double value4 = 0);



    size_t getPointIndex(size_t, size_t);
    std::vector <size_t> getPointCoordinates(size_t index);
    std::vector <size_t> boundary_points();
    Vector Vector_wo_boundary_nodes(Vector &v);

    size_t getN() const;
    size_t getM() const;
    double getHx() const;
    double getHy() const;
};

#endif // GRID_H
