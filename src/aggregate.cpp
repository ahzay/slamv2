//
// Created by user on 7/10/23.
//

#include "aggregate.h"

MatrixXd Aggregate::get_mat() const {
    MatrixXd m;
    m.conservativeResize(data_v.size(), 4);
    for (int i = 0; i < data_v.size(); i++)
        m.row(i) = data_v[i].d;
    return m;
}

void Aggregate::push_back(Data data) {
    data_v.push_back(data);
}

bool Aggregate::check_continuity(Data data) {
    double dist_variation = abs(data.d(2) - data_v.back().d(2));
    double angle_variation = abs(atan2(sin(data.d(3) - data_v.back().d(3)),
                                       cos(data.d(3) - data_v.back().d(3))));
    bool ans = angle_variation > _angle_tolerance || dist_variation > _distance_tolerance;
    if (ans) {
        cout << "Dist variation: " << dist_variation << endl
             << "Angle variation: " << angle_variation << endl
             << "Dist tolerance " << _distance_tolerance << endl
             << "Angle tolerance " << _angle_tolerance << endl;
    }
    return ans;
}

void Aggregate::flush(Data data) {
    data_v.clear();
    data_v.push_back(data);
    _pos = data.p;
}

Aggregate::Aggregate(int angle_tolerance, int distance_tolerance) {
    _angle_tolerance = angle_tolerance;
    _distance_tolerance = distance_tolerance;
}
