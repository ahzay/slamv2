//
// Created by user on 7/12/23.
//

#ifndef MAPPERNODE_IEKF_H
#define MAPPERNODE_IEKF_H

#include "entity.h"
#include "data.h"

class Iekf {
public:
    Iekf(Entity *e);

    void add(const Data &d);

    void update();

    Entity *_e;
    Aggregate a;
};


#endif //MAPPERNODE_IEKF_H
