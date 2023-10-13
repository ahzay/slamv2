//
// Created by user on 7/12/23.
//

#include "iekf.h"

Iekf::Iekf(Entity *e) : _e(e) {}

void Iekf::add(const Data &d) {
    if (_e && _e == d._e)
        a.push_back(d);
}

void Iekf::update() {
    // TODO: make another variable for this ...
    if (a._data_vector.size() >= _e->m->_options.init_npoints/_e->m->_options.init_npoints_multiplier)
        _e->m->ap_ls(*_e, a);
}

