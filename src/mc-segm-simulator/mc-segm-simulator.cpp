//
// Created by user on 10/11/23.
//
#include "simulator.h"
#include <iostream>

int main(int argc, char *argv[]) {
    CmdLineOptions o(argc, argv);
    cout << "to: " << o.data_folder << endl;
    Simulator s(o);
    s.simulate();
}