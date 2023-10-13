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

    ~Entity();

    void print() const;

    double mahalanobis(const Data &data) const;
    pair<double, double> mahalanobis(const Aggregate &aggregate)const;
    bool operator==(const Entity &oe) const;

    void update();

    Entity(const Entity &other);
    Entity& operator=(const Entity &other);


    Model *m;
    VectorXd _p; // parameters
    VectorXd _t; // post processing attributes
    MatrixXd _E; // uncertainty matrix, initialized with DOP
    int updateLeft = 60;
    bool augmented_last_scan=false;
    Aggregate *_a; // test: entity aggregate
};


#endif //MAPPERNODE_ENTITY_H
