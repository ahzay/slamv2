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
    Aggregate(int distance_tolerance, int angle_tolerance);

    MatrixXd get_mat() const;

    void push_back(Data data);

    bool check_continuity(Data data);

    void flush(Data data);

    // private:
    int _angle_tolerance, _distance_tolerance;
    VectorXd _pos;
    vector<Data> data_v;
};


#endif //MAPPERNODE_AGGREGATE_H
