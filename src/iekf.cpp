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
    _e->m->ap_ls(*_e, a);
}

