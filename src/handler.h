//
// Created by user on 7/12/23.
//

#ifndef MAPPERNODE_HANDLER_H
#define MAPPERNODE_HANDLER_H

#include "scan.h"
#include "entitymap.h"
#include "fsm.h"
#include "cmdlineoptions.h"
#include "visualizer.h"
#include <random>

class Handler {
public:
    Handler(const deque<Model *> &models, const CmdLineOptions &options);

    void preprocess_scan(Scan &scan);

    void process_scan();

    State process_measurement(Data &data);

    void f_begin();

    void f_continuous();

    void f_fsm();

    //void f_close_fsms();

    void end(bool at_scan_end);

    // attributes
    const deque<Model *> _models; // permanent
    Data _data;
    unsigned short n = 0;
    int adjusted_init_npoints;
    State state = s_begin;
    EntityMap map;
    Scan ps; // preprocessed scan
    Scan ua; // unassociated aggregate
    /*TODO: implement keeping unassociated and unprocessed data
     * have to figure our how to recheck for continuty (maybe
     * bring all the scans back to one and check continuity)*/
    Scan uua; // unassociated and unused aggregate
    vector<Fsm> fsms;
    vector<Fsm *> fsms_are_done;
    const CmdLineOptions o;
    Visualizer *v;

};


#endif //MAPPERNODE_HANDLER_H
