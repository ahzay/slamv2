//
// Created by user on 10/11/23.
//

#ifndef MC_SEGM_SIMLATOR_CMDLINEOPTIONS_H
#define MC_SEGM_SIMLATOR_CMDLINEOPTIONS_H


#include <string>

using namespace std;

class CmdLineOptions {
public:
    CmdLineOptions(int argc, char *argv[]);

    string data_folder;
    int nell;
    double xmin, xmax, ymin, ymax;
    // an in deg
    double an, npts; // simulated scanner coverage angle and point number per scan
};


#endif //MC_SEGM_SIMLATOR_CMDLINEOPTIONS_H
