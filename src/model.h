//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_MODEL_H
#define MAPPERNODE_MODEL_H

#include "aggregate.h"
#include "tools.h"
#include "cmdlineoptions.h"

class Model {
public:
    explicit Model(const string &file, const CmdLineOptions &options);


    double fs(const Entity &e, const Data &d) const;

    virtual double fss(const VectorXd &p, const VectorXd &pose, const VectorXd &rotated_measurement) const = 0;

    VectorXd fsn(const Entity &e, const Aggregate &a) const;

    virtual Matrix<double, 1, -1> dfs(const Entity &e, const Data &d) const = 0;

    MatrixXd dfsn(const Entity &e, const Aggregate &a) const;

    void dop(Entity &e, const Aggregate &a) const;

    void ap_dop(Entity &e, const Aggregate &a) const; // TODO: this

    virtual VectorXd init(const Aggregate &a) const = 0;

    virtual void ls(Entity &e, const Aggregate &a) const = 0;

    // TODO: generalize this and others so not virtual !!!!
    virtual void ap_ls(Entity &e, const Aggregate &a) const = 0;

    virtual bool safety(Entity &e) const = 0;

    virtual void init_post(Entity &e, const Aggregate &a) const = 0;

    virtual bool associate(const Entity &e, const Data &d) const = 0;

    virtual void augment_post(Entity &e, const Data &d) const = 0;

// private:
    CmdLineOptions _options;
    VectorXd _parameter_mins, _parameter_maxs;
    VectorXi _dfs_parameter_indexes, _dfs_measurement_error_indexes;
    MatrixXd _W_a, _W_s, _Q_a, _Q_s, I;
    double _mahalanobis_strict, _mahalanobis_flex, _dop_sigma,
            _mahalanobis_aug,_mahalanobis_aug_min; // thresholds
    unsigned short _parameter_count, _model_index, _post_attributes_count;
    double _ap_ls_forgetting_factor;
    bool _closed_fs;

};


#endif //MAPPERNODE_MODEL_H
