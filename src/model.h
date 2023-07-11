//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_MODEL_H
#define MAPPERNODE_MODEL_H

#include "aggregate.h"

class Model {
public:
    explicit Model(const string &file);

    virtual VectorXd fs(const Entity &e, const Aggregate &a) = 0;

    virtual MatrixXd dfs(const Entity &e, const Aggregate &a) = 0;

    MatrixXd dop(const Entity &e, const Aggregate &a);

    virtual VectorXd ls(const Aggregate &a) = 0;

    virtual VectorXd ls(Entity &e, const Aggregate &a) = 0;

    virtual bool safety(Entity &e) = 0;

    virtual void init_post(Entity &e, const Aggregate &a) = 0;

    virtual bool associate(const Entity &e, const Data &d) = 0;

    virtual void augment_post(Entity &e, const Data &d) = 0;

// private:
    VectorXd _parameter_mins, _parameter_maxs;
    VectorXi _dfs_parameter_indexes, _dfs_measurement_error_indexes;
    MatrixXd _W_a, _W_s, _Q_a, _Q_s, I;
    double _mahalanobis_strict{}, _mahalanobis_flex{}, _dop_sigma{},
            _mahalanobis_aug{}; // thresholds
    unsigned short _parameter_count{}, _model_index{}, _post_attributes_count{};
    bool _closed_fs{};
};


#endif //MAPPERNODE_MODEL_H
