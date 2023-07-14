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
    Data(Vector2d measurement, Vector3d pose);

    Vector2d get_xy() const;

    Vector2d get_rotated_measurement() const;

    // attributes
    Vector2d _measurement;          // data (measurement)
    Vector3d _pose;          // position at measurement
    Entity *_e = nullptr; // associated entity
};


#endif //MAPPERNODE_DATA_H
