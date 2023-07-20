//
// Created by user on 7/12/23.
//

#include "mappernode.h"

int main(int argc, char *argv[]) {

    CmdLineOptions options(argc, argv);
    deque<Scan> scans(options.scan_number, Scan(options));
    for (int i = 0; i < options.scan_number; i++)
        scans[i].read_scan(i);
    /*auto *ellipse_model1 = new EllipseModel(options.data_folder + "/model_0307_real.txt",
                                            options);
    auto *ellipse_model2 = new EllipseModel(options.data_folder + "/model_0713_real.txt",
                                            options);
    auto *ellipse_model3 = new EllipseModel(options.data_folder + "/model_1317_real.txt",
                                            options);*/
    auto *line_model = new LineModel(options.data_folder + "/model_line_real.txt",
                                     options);
    auto *ellipse_model = new EllipseModel(options.data_folder + "/model_e_real.txt",
                                           options);
    auto e_copy(*ellipse_model);
    deque<Model *> models;
    //models.push_back(ellipse_model1);
    //models.push_back(ellipse_model2);
    //models.push_back(ellipse_model3);
    /*for(double i=0.1;i<1.9;i+=0.1){
        e_copy._parameter_mins[5] = i;
        e_copy._parameter_maxs[5] = i + 0.1;
        models.push_back(new EllipseModel(e_copy));
    }*/
    e_copy._parameter_mins[5] = 0.1;
    e_copy._parameter_maxs[5] = 1.9;
    models.push_back(new EllipseModel(e_copy));
    models.push_back(line_model);
    Handler h(models, options);
    int cnt = 0;
    system("rm *png *gpt");
    for (auto &scan: scans) {
        Visualizer v(cnt++);
        h.v = &v; // IMPORTANT TODO: make this safe
        cout << "*********************************************************" << endl;
        cout << "************************NEW SCAN*************************"
             << cnt - 1 << endl;
        cout << "*********************************************************" << endl;
        h.preprocess_scan(scan);
        h.process_scan(); // also ends it inside
        for (const auto &e: h.map.entities)
            if (e.m->_parameter_count == 6)
                v.add_ellipse(e._p);
            else v.add_segment(e._p, e._t);
        v.save();
        system(("gnuplot " + to_string(cnt - 1) + ".gpt > " +
                to_string(cnt - 1) + ".png").c_str());
    }
    // render plots
    //system("parallel -j 24 gnuplot {} \">\" {.}.png ::: *.gpt");
    return 0;
}