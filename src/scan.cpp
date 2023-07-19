//
// Created by user on 7/11/23.
//

#include "scan.h"
#include <fstream>

bool Scan::check_continuity(Data data) const {
    // TODO: normalize because angles don't mean the same thing at diff distances
    // maybe add eucledian variation test if both other tests pass
    double dist_variation = abs(data._measurement(0) - _data_vector.back()._measurement(0));
    double angle_variation = abs(atan2(sin(data._measurement(1) - _data_vector.back()._measurement(1)),
                                       cos(data._measurement(1) - _data_vector.back()._measurement(1))));
    double eucledian_variation = (data.get_xy() - _data_vector.back().get_xy()).norm();
    bool ans = (angle_variation > _options.angle_tolerance ||
               dist_variation > _options.distance_tolerance) &&
               eucledian_variation > _options.eucledian_tolerance;
    if (ans) {
        cout << "Continuity broken at: " << data.get_xy().transpose() << endl;
        cout << "Dist variation: " << dist_variation << endl
             << "Angle variation: " << angle_variation << endl
             << "Eucledian variation" << eucledian_variation << endl
             << "Dist tolerance " << _options.distance_tolerance << endl
             << "Angle tolerance " << _options.angle_tolerance << endl
             << "Eucledian tolerance" << _options.eucledian_tolerance << endl;
    }
    return ans;
}

Scan::Scan(CmdLineOptions options) : Aggregate(), _options(options) {}

void Scan::read_scan(const int n) {
    char cbuf;
    double dbuf;
    string filename = _options.data_folder + "/combined_" + to_string(n) + ".csv";
    ifstream f(filename);
    if (!f.is_open()) {
        throw std::runtime_error("Error opening the scan file!");
    }
    int rows = 0;
    // first line
    f >> _pose(0) >> cbuf >> _pose(1) >> cbuf >> _pose(2);
    // mult
    _pose({0, 1}) *= _options.env_multiplier;
    // odom x y z q0 q1 q2 q3
    f >> dbuf >> dbuf >> dbuf >> dbuf >> dbuf >> dbuf >> dbuf;
    // ,rows
    f >> cbuf >> rows;
    cout << "pose: " << _pose.transpose() << endl;
    _data_vector.empty();
    // get data
    for (int i = _options.prune; i < rows - _options.prune; i++) {
        Vector<double, 4> bvec;
        f >> bvec(0) >> cbuf >> bvec(1) >> cbuf >> bvec(2) >> cbuf >> bvec(3);
        if (bvec(2) > 0 && bvec(2) < _options.max_scan_distance) { // restriction on min/max distance
            push_back(Data(bvec({2, 3}), _pose, _options.data_longevity));
        }
    }
    // mult and error
    for (auto &d: _data_vector) {
        double error = (-0.01 + (rand() / (RAND_MAX / 0.02)));
        error = 0;
        d._measurement(0) += error;
        d._measurement(0) *= _options.env_multiplier;
    }
    /*
    //  adjust data to emulate null rotation
    data.col(3).array() += ori;
    data.col(0) = data.col(2).array() * data.col(3).array().cos() + loc(0);
    data.col(1) = data.col(2).array() * data.col(3).array().sin() + loc(1);*/
}
