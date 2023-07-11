//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_ENTITY_H
#define MAPPERNODE_ENTITY_H

#include "model.h"
class Entity {
public:
    Entity(Model* model, VectorXd p, MatrixXd E);
    double mahalanobis(Data& data);
    Model* m;
    VectorXd _p; // parameters
    VectorXd _t; // post processing attributes
    MatrixXd _E; // uncertainty matrix, initialized with DOP
};


#endif //MAPPERNODE_ENTITY_H
