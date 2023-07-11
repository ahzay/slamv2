//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_DATA_H
#define MAPPERNODE_DATA_H
#include "Eigen/Core"
using namespace Eigen;
class Entity;
class Data {
public:
    Data() = default;
    Data(VectorXd data, VectorXd pose);
    // attributes
    VectorXd d;          // data (measurement)
    VectorXd p;          // pose at measurement
    Entity *e = nullptr; // associated entity
};


#endif //MAPPERNODE_DATA_H
