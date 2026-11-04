#ifndef BC_H
#define BC_H

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Core"
#include "settings.h"

class Bc
{
private:
    enum BCnames : std::size_t{
        BCx,
        BCy,
        Fx,
        Fy,
        countBCnames
    };

    std::array<std::vector<int>,countBCnames> bcs {};

public:
    Bc() = default;

    void extractBCSalome();
    void extractBCAbaqus();

    void setBc(Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& disp
        ,Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& force) const;
};

#endif