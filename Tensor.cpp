#include "Tensor.h"


Tensor::Tensor(): N(10), M(10), hx(1/9), hy(1/9), k_r(1), k_t(1)
{

}

Tensor::Tensor(size_t n, size_t m, double h_x, double h_y, double kr, double kt): N(n), M(m), hx(h_x), hy(h_y), k_r(kr), k_t(kt)
{

}

double Tensor::K_11(size_t i, size_t j)
{
    double x = i * hx;
    double y = j * hy;
    return (k_r * x * x + k_t * y * y)/(x * x + y * y);
}

double Tensor::K_12(size_t i, size_t j)
{
    double x = i * hx;
    double y = j * hy;
    return  0.0;
    return (k_r - k_t) * x * y/(x * x + y * y);
}

double Tensor::K_21(size_t i, size_t j)
{
    double x = i * hx;
    double y = j * hy;
    return  0.0;
    return (k_r - k_t) * x * y/(x * x + y * y);
}

double Tensor::K_22(size_t i, size_t j)
{
    double x = i * hx;
    double y = j * hy;
    return (k_r * y * y + k_t * x * x)/(x * x + y * y);
}

double Tensor::getValue(size_t i, size_t j, int num_elem)
{
    if (num_elem == 1)
        return K_11(i, j);
    if (num_elem == 2)
        return K_12(i, j);
    if (num_elem == 3)
        return K_21(i, j);
    else // num_elem = 4
        return K_22(i, j);
}
