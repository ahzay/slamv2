//
// Created by user on 10/11/23.
//

#include "cmdlineoptions.h"
#include <iostream>

CmdLineOptions::CmdLineOptions(int argc, char **argv) {
    if (argc < 8)
        throw std::runtime_error("Not enough input arguments.");
    xmin = stod(argv[1]);
    xmax = stod(argv[2]);
    ymin = stod(argv[3]);
    ymax = stod(argv[4]);
    nell = stoi(argv[5]);
    an = stod(argv[6]);
    npts = stod(argv[7]);
    data_folder = argv[8];
}
