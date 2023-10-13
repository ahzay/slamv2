//
// Created by user on 7/11/23.
//

#ifndef MAPPERNODE_SCAN_H
#define MAPPERNODE_SCAN_H

#include "aggregate.h"
#include "cmdlineoptions.h"
#include "tools.h"
class Scan : public Aggregate {
public:
    explicit Scan(CmdLineOptions options);

    void read_scan(int n);

    bool check_continuity(Data data) const;

    CmdLineOptions _options;
    //Vector3d _pose; // x,y, rot
};


#endif //MAPPERNODE_SCAN_H
