#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>

class Tensor
{
    size_t N, M;
    double hx, hy;
    double k_r, k_t;
public:
    Tensor();
    Tensor(size_t, size_t, double, double, double, double);

    double K_11(size_t, size_t); // first tensor element (0, 0)
    double K_12(size_t, size_t);
    double K_21(size_t, size_t);
    double K_22(size_t, size_t);

    double getValue(size_t, size_t, int);
};

#endif // TENSOR_H
