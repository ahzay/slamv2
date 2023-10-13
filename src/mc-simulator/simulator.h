//
// Created by user on 10/11/23.
//

#ifndef MAPPERNODE_SIMULATOR_H
#define MAPPERNODE_SIMULATOR_H

#include "../mappernode/tools.h"
#include "cmdlineoptions.h"
#include <random>

using namespace Eigen;

class Simulator {
public:
    Simulator(CmdLineOptions options);

    double get_rand(double min, double max);

    void gen_params();

    void gen_location();

    void simulate();

    CmdLineOptions o;
    mt19937 generator;
    Vector<double, 6> p; // parameters
    Vector2d l; // robot location
    vector<Vector2d> ms; // measurements
};


#endif //MAPPERNODE_SIMULATOR_H
