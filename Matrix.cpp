#include "Matrix.h"

Matrix::Matrix() : N(10), M(10)
{
    setSize();
}

Matrix::Matrix(size_t rows, size_t cols) : N(rows), M(cols)
{
    setSize();
}

size_t Matrix::getN() const
{
    return N;
}

size_t Matrix::getM() const
{
    return M;
}



void Matrix::setSize()
{
    field.resize(M);
    for (size_t i = 0; i < M; i++)
        field[i].resize(N);
}

void Matrix::setValue(size_t x, size_t y, double value)
{
    field[x][y] = value;
}

double Matrix::getValue(size_t x, size_t y) const
{
    return field[x][y];
}


double Matrix::max()
{
    double c = 0.0;

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < M; j++)
            if (abs(this->getValue(i,j)) > c)
                c = abs(this->getValue(i, j));
    return c;
}

void Matrix::swapRows(size_t k, size_t l)
{
    double temp;

    for (size_t j = 0; j < M; j++) {
       temp = field[k][j];
       field[k][j] = field[l][j];
       field[l][j] = temp;
    }
}

std::ostream & operator<<(std::ostream & os, const Matrix & m)
{
    for (size_t i = 0; i < m.N; i++) {
        for (size_t j = 0; j < m.M; j++)
            os << m.field[i][j] << " ";
    os << std::endl;
    }
    return os;
}

void Matrix::fillNonZeroElementsPosition()
{
    nonZeroElemntsPosition.resize(N);
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < M; j++)
            if(fabs(getValue(i,j)) > eps)
                nonZeroElemntsPosition[i].push_back(j);
        }
}

void Matrix::saveToFile(std::string filename) const
{
    std::ofstream fout(filename);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++)
            fout << field[i][j] << ' ';
        fout << std::endl;
    }

    fout.close();
}


void Matrix::copyFromData(const Matrix &m)
{
    if (this->getM() != m.getM()) {
        std::cout << "Function copyFromData: sizes aren't equale" << std::endl;
        return;
    }

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < M; j++)
            setValue(i, j, m.getValue(i,j));
}

void trans(Matrix &D1, Matrix &D2)
{
    for (size_t i = 0; i < D1.getN(); i++) {
        for (size_t j = 0; j < D1.getM(); j++) {
            D2.setValue(i, j, D1.getValue(j, i));
        }
    }
}

void sqr(Matrix &D1)
{
    double c = 0.0;
    for(size_t i = 0; i < D1.getN(); i++)
        for (size_t j = 0; j < D1.getM(); j++)
        {
            c = D1.getValue(i, j) * D1.getValue(i, j);
            D1.setValue(i, j, c);

        }
}
