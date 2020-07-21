#include "Setting.h"


Setting::Setting(size_t N_, size_t M_): N(N_), M(M_)
{

}
/*
void Setting::fileInitialization(std::string filename)
{
    std::ofstream fout(filename);

    setGrid(fout);
    setTensor(fout);
    setMatrix(fout);
    setRightPart(fout);

    fout.close();
}
*/
void Setting::setGrid(std::ofstream &fout)
{
    std::vector <double> steps(2);

    steps = setStep();
    fout << N << ' ' << M << ' ' << steps[0] << ' ' << steps[1] << std::endl;
}


std::vector<double> Setting::setStep()
{
    std::vector <double> steps(2);

    double hx, hy;

    hx = 1.0 / N;
    hy = 1.0 / M;

    steps[0] = hx;
    steps[1] = hy;

    return steps;
}
