//
// Created by user on 7/10/23.
//

#include "data.h"

Data::Data(Vector2d measurement, Vector3d pose,
           int longevity) : _measurement(measurement),
                            _pose(pose), life(longevity) {}

Vector2d Data::get_rotated_measurement() const {
    return Vector2d(_measurement(0), _measurement(1) + _pose(2));
}

void Data::set_xy(Vector2d xy) { // make sure pose is set accurately first
    auto dx = xy(0) - _pose(0);
    auto dy = xy(1) - _pose(1);
    _measurement(0) = sqrt(dx * dx + dy * dy);
    auto an = atan2(dy, dx);
    _measurement(1) = an - _pose(2);
}

Vector2d Data::get_xy() const {
    auto dan = get_rotated_measurement();
    return {dan(0) * cos(dan(1)) + _pose(0),
            dan(0) * sin(dan(1)) + _pose(1)};
}

void Data::change_referential(Vector3d new_pose) {
    if (new_pose != _pose) {
        auto dxdy = get_xy() - new_pose({0, 1});
        _pose = new_pose;
        auto rotated_angle = atan2(dxdy(1), dxdy(0));
        auto d = sqrt(dxdy(1) * dxdy(1) + dxdy(0) * dxdy(0));
        _measurement(0) = d;
        _measurement(1) = rotated_angle - _pose(2);
    }
}

void Data::normalize() {
    double angle = std::fmod(_measurement(1), 2.0 * M_PI);
    if (angle < 0)
        angle += 2.0 * M_PI;
    _measurement(1) = angle;
}

double Data::distance_to(const Data &o) const {
    return (get_xy() - o.get_xy()).norm();
}

