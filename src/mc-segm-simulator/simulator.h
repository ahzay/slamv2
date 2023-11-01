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

    void gen_all();

    void simulate();

    CmdLineOptions o;
    mt19937 generator;
    vector<Vector<double, 6>> ps; // parameters
    Vector2d l; // robot location
    vector<Vector2d> ms; // measurements

private:
    double get_rand(double min, double max);

    Vector<double, 6> gen_params();

    bool detect_overlap(Vector<double, 6> p1, Vector<double, 6> p2);

    void gen_all_params();

    void gen_location();

    double calc_e(Vector<double,6> p, double x, double y);
};


#endif //MAPPERNODE_SIMULATOR_H
