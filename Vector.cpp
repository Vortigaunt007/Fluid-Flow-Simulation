#include "Vector.h"

Vector::Vector(): M(10)
{
    setSize();
}

Vector::Vector(size_t cols): M(cols)
{
    setSize();
}

void Vector::setSize()
{
    str.resize(M);
}

size_t Vector::getSize() const
{
    return M;
}

double &Vector::operator[](size_t i)
{
      return str[i];
}

double Vector::operator[](size_t i) const
{

    return str[i];
}

Vector &Vector::operator =(const Vector &v)
{
    for (size_t i = 0; i < M; i++)
        str[i] = v[i];
    return *this;
}

void Vector::saveToFile(size_t num_string, size_t num_column, std::string filename) const
{
   // std::ofstream fout(filename);

    Matrix B_grid(num_string, num_column);

    int num_s, num_c;
    size_t count = 0;
    num_s = num_string - 1;
    num_c = 0;

    for (size_t i = 0; i < this->getSize(); i++) {

        //std::cout << "Count = " << count << " " << num_s << " " << num_c << " "  << std::endl;
        count++;
        B_grid.setValue(num_s, num_c, str[i]);
        num_s--;
        if (num_s < 0)
            num_s = num_string - 1;
        if (count != 1 && (count % num_string ) == 0)
            num_c++;
    }

    B_grid.saveToFile(filename);
}

std::ostream & operator<<(std::ostream &os, const Vector &s)
{
    for (size_t i = 0; i < s.M; i++)
            os << s.str[i] << " ";
    os << std::endl;

    return os;
}


double Vector::normL2()
{
    return sqrt(dotProduct(*this, *this));
}


bool checkSizes(const Vector &v1, const Vector &v2)
{
    return !(v1.getSize() != v2.getSize());
}


double dotProduct(const Vector & v1,const Vector & v2)
{
    double result = 0.0;

    if (!checkSizes(v1, v2)) {
        std::cout << "dotProduct: sizes aren't equale" << std::endl;
        return result;
    }

    for (size_t i = 0; i < v1.getSize(); i++)
            result += v1[i] * v2[i];

    return result;
}

void summAlpha(const Vector &v1, double alpha, const Vector &v2, Vector &out)
{
    if(!checkSizes(v1, v2)) {
        std::cout << "Function summAlpha: sizes aren't equale" << std::endl;
        return;
    }

    double value = 0.0;
    for(size_t i = 0; i < v1.getSize(); i++) {
        value = v1[i] + alpha * v2[i];
        out[i] = value;
    }
}

void mult_Matrix_Vector(Matrix &D, const Vector &S, Vector &out)
{

    if(D.getM() != S.getSize()){
        std::cout << "Function mult_Matrix_Vector: sizes aren't equale" << std::endl;
        return;
    }

    for(size_t i = 0; i < D.getN(); i++) {
            double res = 0.0;
            for(size_t k : D.nonZeroElemntsPosition[i])
                res += D.getValue(i, k) * S[k];
            out[i] = res;
    }
}


void sqr(Vector &v)
{
    double c = 0.0;

    for(size_t i = 0; i < v.getSize(); i++) {
            c = v[i] * v[i];
            v[i] = c;
    }
}
