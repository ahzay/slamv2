//
// Created by user on 7/12/23.
//

#ifndef MAPPERNODE_FSM_H
#define MAPPERNODE_FSM_H

#include "model.h"
#include "aggregate.h"
#include "entity.h"
#include "ekf.h"

class Fsm {
public:
    Fsm(Model *model, const Aggregate &a);

    void f_flexible();

    void f_strict();

    void f_least_squares(const bool is_scan_end);

    State process_measurement(Data &data);

    //bool operator<(const Fsm &other_fsm) const;

    Data _data;
    Aggregate _a;
    Entity e;
    State state = s_fsm_flexible;
    Ekf ekf;
};


#endif //MAPPERNODE_FSM_H
