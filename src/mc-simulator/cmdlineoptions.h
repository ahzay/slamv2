//
// Created by user on 10/11/23.
//

#ifndef MAPPERNODE_CMDLINEOPTIONS_H
#define MAPPERNODE_CMDLINEOPTIONS_H


#include <string>

using namespace std;

class CmdLineOptions {
public:
    CmdLineOptions(int argc, char *argv[]);

    string data_folder;
    double xmin, xmax, ymin, ymax;
    // an in deg
    double an, npts; // simulated scanner coverage angle and point number per scan
};


#endif //MAPPERNODE_CMDLINEOPTIONS_H
