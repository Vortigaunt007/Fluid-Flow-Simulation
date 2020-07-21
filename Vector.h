#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>
#include <math.h>
#include "Matrix.h"

class Vector
{
private:
    std::vector<double> str;
    size_t M;

    void setSize();
public:
    Vector();
    Vector(size_t);

    size_t getSize() const;

    double& operator[](size_t i);
    double operator[](size_t i) const;

    Vector &operator =(const Vector &v);
    friend std::ostream& operator<<(std::ostream&, const Vector&);

    void saveToFile(size_t num_string, size_t num_column, std::string filename) const;

    double normL2();
};
void sqr(Vector &v);

bool checkSizes(const Vector&, const Vector&);
double dotProduct(const Vector&,const Vector&);
void summAlpha(const Vector&, double, const Vector&, Vector&);
void mult_Matrix_Vector(Matrix&, const Vector&, Vector&);

#endif // VECTOR_H
