#ifndef SETTING_H
#define SETTING_H

#include <iostream>
#include <fstream>
#include <vector>

class Setting
{
    size_t N, M;
public:
    Setting(size_t N_, size_t M_);
    void fileInitialization(std::string filename);
    void setGrid(std::ofstream &fout);
    void setTensor(std::ofstream &fout);
    void setMatrix(std::ofstream &fout);
    void setRightPart(std::ofstream &fout);

    std::vector <double> setStep();
};

#endif // SETTING_H
