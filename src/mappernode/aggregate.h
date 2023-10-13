//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_AGGREGATE_H
#define MAPPERNODE_AGGREGATE_H

#include "data.h"
#include <vector>
#include <iostream>

using namespace std;

class Aggregate {
public:
    Aggregate();

    //Aggregate(const MatrixX2d &mat, const Vector3d &pose);

    MatrixX2d get_measurement_mat() const;

    MatrixX2d get_rotated_measurement_mat() const;

    MatrixX2d get_xy_mat() const;

    void change_referential(Vector3d new_pose);

    void self_sort();

    void reorder();

    void push_back(const Data &data);

    void push_back(const Aggregate &a);

    void clear();

    void flush(const Data &data);

    void flush();

    vector<Data> _data_vector; // non rotated data
    Vector3d _pose;
};


#endif //MAPPERNODE_AGGREGATE_H
