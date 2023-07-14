//
// Created by user on 7/10/23.
//

#include "data.h"

Data::Data(Vector2d measurement, Vector3d pose) : _measurement(measurement),
                                                  _pose(pose) {}

Vector2d Data::get_rotated_measurement() const {
    return Vector2d(_measurement(0), _measurement(1) + _pose(2));
}

Vector2d Data::get_xy() const {
    auto dan = get_rotated_measurement();
    return Vector2d(dan(0) * cos(dan(1)) + _pose(0),
                    dan(0) * sin(dan(1)) + _pose(1));
}
