//
// Created by user on 7/12/23.
//

#ifndef MAPPERNODE_ENTITYMAP_H
#define MAPPERNODE_ENTITYMAP_H

#include "entity.h"
#include "iekf.h"
#include <deque>

class EntityMap {
public:
    EntityMap();

    bool associate_data(Data &d);

    void add_augment(const Data &data);

    void run_augment();

    void add_entity(const Entity &e, const Aggregate &a);


    deque<Entity> entities;
    deque<Iekf> iekfs;
};


#endif //MAPPERNODE_ENTITYMAP_H
