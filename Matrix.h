#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>

class Matrix
{
private:
    std::vector<std::vector<double>> field;
    size_t N, M; // string, column
    double eps = 0.000001; // fillNonZeroElementsPosition

   void setSize();
public:
    std::vector<std::vector<size_t>> nonZeroElemntsPosition;

    Matrix();
    Matrix(size_t, size_t);

    void copyFromData(const Matrix &);

    void setValue(size_t, size_t, double);
    double getValue(size_t, size_t) const;

    size_t getN() const;
    size_t getM() const;

    void fillNonZeroElementsPosition();


    void saveToFile(std::string filename) const;

    friend std::ostream& operator<<(std::ostream&, const Matrix&);


    double max();
    void swapRows(size_t, size_t);
};

void trans(Matrix &D1, Matrix &D2);
void sqr(Matrix &D1);


#endif // MATRIX_H
