//
// Created by user on 7/11/23.
//

#ifndef MAPPERNODE_CMDLINEOPTIONS_H
#define MAPPERNODE_CMDLINEOPTIONS_H

#include <string>

using namespace std;

class CmdLineOptions {
public:
    CmdLineOptions(int argc, char *argv[]);

    string data_folder;
    int scan_number, stride, prune, init_npoints, data_longevity,npoints_mult;
    double angle_tolerance, distance_tolerance, env_multiplier, max_scan_distance,eucledian_tolerance;
};


#endif //MAPPERNODE_CMDLINEOPTIONS_H
