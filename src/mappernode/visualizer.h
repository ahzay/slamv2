//
// Created by user on 7/13/23.
//

#ifndef MAPPERNODE_VISUALIZER_H
#define MAPPERNODE_VISUALIZER_H

#include <Eigen/Core>
#include <fstream>
#include "aggregate.h"

using namespace std;
using namespace Eigen;

class Visualizer {
public:
    Visualizer(int cnt);

    void add_points(Eigen::Matrix<double, Eigen::Dynamic, 2> pts, const string &color);

    void add_data(const Data &d, const string &color);

    void add_aggregate(const Aggregate &a, const string &color);

    void add_ellipse(Eigen::VectorXd p, const std::string &color);

    void add_segment(Eigen::VectorXd p, Eigen::VectorXd t);

    void save();

    ofstream of;
};


#endif //MAPPERNODE_VISUALIZER_H
