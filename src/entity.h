//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_ENTITY_H
#define MAPPERNODE_ENTITY_H

#include "model.h"

class Entity {
public:
    Entity(Model *model, const VectorXd &p, const MatrixXd &E);

    Entity(Model *model, const Aggregate &a);

    void print() const;

    double mahalanobis(const Data &data) const;

    bool operator==(const Entity &oe) const;

    Model *m;
    VectorXd _p; // parameters
    VectorXd _t; // post processing attributes
    MatrixXd _E; // uncertainty matrix, initialized with DOP
};


#endif //MAPPERNODE_ENTITY_H
