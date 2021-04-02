#include "Grid.h"


size_t Grid::getN() const
{
    return N;
}

size_t Grid::getM() const
{
    return M;
}

double Grid::getHx() const
{
    return hx;
}

double Grid::getHy() const
{
    return hy;
}

Grid::Grid(): N(10), M(10), hx(1/9), hy(1/9)
{

}

Grid::Grid(size_t n, size_t m): N(n), M(m)//, hx(h_x), hy(h_y)
{
    setGridStep();
}

void Grid::setGridStep()
{
    hx = 1.0 / (N - 1);
    hy = 1.0 / (M - 1);
}

Grid::typePoint Grid::type(size_t i, size_t j)
{
    if ((i == 0 && j == 0) ||
            (i == 0 && j == N - 1) ||
            (i == M - 1 && j == 0) ||
            (i == M - 1 && j == N - 1))
            return ANGLE;
        if (i == 0) return LEFT;
        if (i == M - 1) return RIGHT;
        if (j == 0) return BOTTOM;
        if (j == N - 1) return TOP;
        return CENTER;
}

Matrix Grid::initialization_A_test()
{
    Grid A_Grid(N, M);
    Matrix A(M * N, M * N);
    int line_index = 0;

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            line_index = A_Grid.getPointIndex(i, j);
            A.setValue(line_index, A_Grid.getPointIndex(i, j), line_index + 1);
        }
    }

 //   std::cout << A << std::endl;

    return A;
}

Vector Grid::initialization_b_test()
{
     Vector b(N * M);

     for (size_t i = 0; i < M * N; i++) {
         b[i] = i + 1.0;
     }

     //std::cout << b << std::endl;

     return b;
}


double Grid::coef_a_x(size_t i, size_t j, Tensor &K_Tensor)
{
    return 2 * K_Tensor.K_11(i - 1, j) * K_Tensor.K_11(i, j)
            / (K_Tensor.K_11(i - 1, j) + K_Tensor.K_11(i, j));
}


double Grid::coef_a_y(size_t i, size_t j, Tensor &K_Tensor)
{
    return 2 * K_Tensor.K_22(i, j - 1) * K_Tensor.K_22(i, j)
            / (K_Tensor.K_22(i, j - 1) + K_Tensor.K_22(i, j));
}


Matrix Grid::initialization_A(Tensor &K_Tensor)
{
    size_t line_index = 0;
    typePoint node_type;
    Grid A_Grid(N, M);
    Matrix A(M * N, M * N);

    // first order of approximation
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            line_index = A_Grid.getPointIndex(i, j);
            node_type = A_Grid.type(i, j);
            if (node_type == ANGLE) {
                //std::cout << "ANGLE" << std::endl;
                A.setValue(line_index, i * N + j, 1.0);
            } else if (node_type == LEFT) {
                //std::cout << "LEFT"  << std::endl;

               /* // first order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(1, j), 1.0 / hx * K_Tensor.K_11(1,j));
                A.setValue(line_index, A_Grid.getPointIndex(0, j), -1.0 / hx * K_Tensor.K_11(0,j));
               */
                // second order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(0, j), -1.0 / hx + (K_Tensor.K_11(1, j) - K_Tensor.K_11(0, j))/(2 * K_Tensor.K_11(0, j) * hx));
                A.setValue(line_index, A_Grid.getPointIndex(1, j), 1.0 / hx + hx / (hy * hy) * (K_Tensor.K_22(0, j) / K_Tensor.K_11(0, j))
                                                                       - (K_Tensor.K_11(1, j) - K_Tensor.K_11(0, j))/(2 * K_Tensor.K_11(0, j) * hx)
                                                                       + hx * (K_Tensor.K_22(1, j + 1) - K_Tensor.K_22(1, j)) / (2 * hy * hy * K_Tensor.K_11(0, j)));
                A.setValue(line_index, A_Grid.getPointIndex(1, j + 1), - hx / (2 * hy * hy) * (K_Tensor.K_22(0, j) / K_Tensor.K_11(0, j)) - hx * (K_Tensor.K_22(1, j + 1) - K_Tensor.K_22(1, j)) / (2 * hy * hy * K_Tensor.K_11(0, j)));
                A.setValue(line_index, A_Grid.getPointIndex(1, j - 1), - hx / (2 * hy * hy) * (K_Tensor.K_22(0, j) / K_Tensor.K_11(0, j)));

            } else if (node_type == RIGHT) {
                //std::cout << "RIGHT" << std::endl;
              /*  // first order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(N - 1, j), 1.0 / hx * K_Tensor.K_11(N - 1,j));
                A.setValue(line_index, A_Grid.getPointIndex(N - 2, j), -1.0 / hx * K_Tensor.K_11(N - 2,j));                
               */
                // second order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(N - 1, j), -1.0 / hx + (K_Tensor.K_11(N - 2, j) - K_Tensor.K_11(N - 1, j))/(2 * K_Tensor.K_11(N - 1, j) * hx));
                A.setValue(line_index, A_Grid.getPointIndex(N - 2, j), 1.0 / hx + hx / (hy * hy) * (K_Tensor.K_22(N - 1, j) / K_Tensor.K_11(N - 1, j))
                                                                       - (K_Tensor.K_11(N - 2, j) - K_Tensor.K_11(N - 1, j))/(2 * K_Tensor.K_11(N - 1, j) * hx)
                                                                       + hx * (K_Tensor.K_22(N - 2, j + 1) - K_Tensor.K_22(N - 2, j)) / (2 * hy * hy * K_Tensor.K_11(N - 1, j)));
                A.setValue(line_index, A_Grid.getPointIndex(N - 2, j + 1), - hx / (2 * hy * hy) * (K_Tensor.K_22(N - 1, j) / K_Tensor.K_11(N - 1, j)) - hx * (K_Tensor.K_22(N - 2, j + 1) - K_Tensor.K_22(N - 2, j)) / (2 * hy * hy * K_Tensor.K_11(N - 1, j)));
                A.setValue(line_index, A_Grid.getPointIndex(N - 2, j - 1), - hx / (2 * hy * hy) * (K_Tensor.K_22(N - 1, j) / K_Tensor.K_11(N - 1, j)));

            } else if (node_type == BOTTOM) {
                //std::cout << "BOTTOM" << std::endl;
                // first order of approximation
               /* A.setValue(line_index, A_Grid.getPointIndex(i, 1), 1.0 / hy * K_Tensor.K_22(i,1));
                A.setValue(line_index, A_Grid.getPointIndex(i, 0), -1.0 / hy * K_Tensor.K_22(i,0));
                */
                // second order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(i, 0), -1.0 / hy + (K_Tensor.K_22(i, 1) - K_Tensor.K_22(i, 0))/(2 * K_Tensor.K_22(i, 0) * hy));
                A.setValue(line_index, A_Grid.getPointIndex(i, 1), 1.0 / hy + hy / (hx * hx) * (K_Tensor.K_11(i, 0) / K_Tensor.K_22(i, 0))
                                                                       + hy * (K_Tensor.K_11(i + 1, 1) - K_Tensor.K_11(i, 1))/(2 * hx * hx * K_Tensor.K_22(i, 0))
                                                                       - (K_Tensor.K_22(i, 1) - K_Tensor.K_22(i, 0)) / (2 * hy * K_Tensor.K_22(i, 0)));
                A.setValue(line_index, A_Grid.getPointIndex(i + 1, 1), - hy / (2 * hx * hx) * (K_Tensor.K_11(i, 0) / K_Tensor.K_22(i, 0)) - hy * (K_Tensor.K_11(i + 1, 1) - K_Tensor.K_11(i, 1)) / (2 * hx * hx * K_Tensor.K_22(i, 0)));
                A.setValue(line_index, A_Grid.getPointIndex(i - 1, 1), - hy / (2 * hx * hx) * (K_Tensor.K_11(i, 0) / K_Tensor.K_22(i, 0)));

            } else if (node_type == TOP) {
                //std::cout << "TOP" << std::endl;
                // first order of approximation
              /*  A.setValue(line_index, A_Grid.getPointIndex(i, N - 1), 1.0 / hy * K_Tensor.K_22(i,N - 1));
                A.setValue(line_index, A_Grid.getPointIndex(i, N - 2), -1.0 / hy * K_Tensor.K_22(i,N - 2));
                */
                // second order of approximation
                A.setValue(line_index, A_Grid.getPointIndex(i, N - 1), -1.0 / hy + (K_Tensor.K_22(i, N - 2) - K_Tensor.K_22(i, N - 1))/(2 * K_Tensor.K_22(i, N - 1) * hy));
                A.setValue(line_index, A_Grid.getPointIndex(i, N - 2), 1.0 / hy + hy / (hx * hx) * (K_Tensor.K_11(i, N - 1) / K_Tensor.K_22(i, N - 1))
                                                                       + hy * (K_Tensor.K_11(i + 1, N - 2) - K_Tensor.K_11(i, N - 2))/(2 * hx * hx * K_Tensor.K_22(i, N - 1))
                                                                       - (K_Tensor.K_22(i, N - 2) - K_Tensor.K_22(i, N - 1)) / (2 * hy * K_Tensor.K_22(i, N - 1)));
                A.setValue(line_index, A_Grid.getPointIndex(i + 1, N - 2), - hy / (2 * hx * hx) * (K_Tensor.K_11(i, N - 1) / K_Tensor.K_22(i, N - 1)) - hy * (K_Tensor.K_11(i + 1, N - 2) - K_Tensor.K_11(i, N - 2)) / (2 * hx * hx * K_Tensor.K_22(i, N - 1)));
                A.setValue(line_index, A_Grid.getPointIndex(i - 1, N - 2), - hy / (2 * hx * hx) * (K_Tensor.K_11(i, N - 1) / K_Tensor.K_22(i, N - 1)));

            } else if (node_type == CENTER) {
                //std::cout << "CENTER" << std::endl;

                double h_y, h_x;
                h_y = 1.0 / hy / hy;
                h_x = 1.0 / hx / hx;

                A.setValue(line_index, A_Grid.getPointIndex(i, j - 1), h_y * coef_a_y(i, j,  K_Tensor));
                A.setValue(line_index, A_Grid.getPointIndex(i - 1, j), h_x * coef_a_x(i, j,  K_Tensor));
                A.setValue(line_index, A_Grid.getPointIndex(i, j), - h_x * coef_a_x(i + 1, j,  K_Tensor) - h_x * coef_a_x(i, j,  K_Tensor)
                                                              - h_y * coef_a_y(i, j + 1,  K_Tensor) - h_y * coef_a_y(i, j,  K_Tensor));
                A.setValue(line_index, A_Grid.getPointIndex(i + 1, j), h_x * coef_a_x(i + 1, j,  K_Tensor));
                A.setValue(line_index, A_Grid.getPointIndex(i, j + 1), h_y * coef_a_y(i, j + 1,  K_Tensor));
            }
        }
    }

   // std::cout << A << std::endl;

    return A;
}

Vector Grid::initialization_b(Tensor &K, Matrix &A, size_t index_dirichlet, double value_dirichlet,
                              size_t l1_start, size_t l1_end, double value1,
                              size_t l2_start, size_t l2_end, double value2,
                              size_t l3_start, size_t l3_end, double value3,
                              size_t l4_start, size_t l4_end, double value4)
{

    Vector b(N * M);    
    Grid A_Grid(N, M);
    double value_temp;

    size_t index = 0;

    for (size_t i = 0; i < b.getSize(); i++)
        b[i] = 0;

    //условия Неймана
    for (size_t i = l1_start; i < l1_end; i++) {
        value_temp = (l1_end - l1_start) * getHx() * K.K_11(0, i) * value1;
        index = A_Grid.getPointIndex(0, i);
        b[index] = value_temp;
    }

    for (size_t i = l2_start; i < l2_end; i++) {
        value_temp = (l2_end - l2_start) * getHx() * K.K_11(M - 1, i) * value2;
        index = A_Grid.getPointIndex(M - 1, i);
        b[index] = value_temp;
    }

    for (size_t i = l3_start; i < l3_end; i++) {
        value_temp = value3 * (l3_end - l3_start) * getHx() / K.K_22(i, 0);
        index = A_Grid.getPointIndex(i ,0);
        b[index] = value_temp;
    }

    for (size_t i = l4_start; i < l4_end; i++) {
        value_temp = value4 * (l4_end - l4_start) * getHx() / K.K_22(i, N - 1) ;
        index = A_Grid.getPointIndex(i, N - 1);
        b[index] = value_temp;
    }

    // условие Дирихле
    for (size_t j = 0; j < A.getM(); j++)
        A.setValue(index_dirichlet, j, 0.0);
    A.setValue(index_dirichlet, index_dirichlet, 1.0);
    b[index_dirichlet] = value_dirichlet;

    return b;
}

size_t Grid::getPointIndex(size_t i, size_t j)
{
    return i * N + j;
}

std::vector <size_t> Grid::getPointCoordinates(size_t index)
{
    std::vector <size_t> IJ(2);

    IJ[0] = index / N;
    IJ[1] = index % N;

    return IJ;
}

std::vector <size_t> Grid::boundary_points()
{
    std::vector <size_t> indices;
    size_t index;

    for (size_t i = 0; i < getM(); i++)
        for (size_t j = 0; j < getN(); j++) {
            typePoint index_type = type(i, j);
            if (index_type == LEFT || index_type == RIGHT || index_type == BOTTOM || index_type == TOP || index_type == ANGLE) {
                index = getPointIndex(i, j);
                indices.push_back(index);
            }
        }

    return indices;
}



Vector Grid::Vector_wo_boundary_nodes(Vector &v)
{
    std::vector <size_t> boundary_points_indices = boundary_points();
    Vector res(v.getSize() - boundary_points_indices.size());
    bool flag = false;
    int count = 0;

    for (size_t i = 0; i < getM(); i++) {
        for (size_t j = 0; j < getN(); j++) {
            flag = false;
            size_t index = getPointIndex(i, j);
            for (size_t k = 0; k < boundary_points_indices.size(); k++)
                if (index == boundary_points_indices[k]) { // узел граничный
                    flag = true;
                    break;
            }
            if (!flag) // узел не граничный
                res[count++] = v[index];
        }
    }

    return res;
}
