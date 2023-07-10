#include "handler.hpp"
#include <iostream>

using namespace std;

/*void simulate() {
  cout << "Define environment and path then simulate" << endl;
  vector<Object> os = {
      {{6, 4, M_PI / 4, 2, 2, 1.2}, 0},             //
      {{M_PI / 2, 6, M_PI / 4, 3 * M_PI / 4}, 1},   //
      {{3 * M_PI / 4, 8.5, 3 * M_PI / 4, M_PI}, 1}, //
      {{10, 10, M_PI / 4, 2, 2, 1.9}, 0},
      //{{M_PI, 12.10, M_PI, 5 * M_PI / 4}, 1},       //
      //{{0, 0, M_PI / 4, 2, 2, 0.7}, 0},
      //{{M_PI / 4, 6, 0, M_PI / 2}, 1}, // line
  };
  for (auto &o : os)
    v.add_object(o, "r-");
  vector<vecext<double>> ls;
  ls.push_back({0, 0});
  ls.push_back({1, 0});
  ls.push_back({2, 0});
  ls.push_back({3, 0});
  ls.push_back({4, 0});
  ls.push_back({5, 0});
  ls.push_back({6, 0});
  ls.push_back({7, 0});
  ls.push_back({8, 0});
  ls.push_back({9, 0});
  ls.push_back({10, 0});
  ls.push_back({10, 1});
  ls.push_back({10, 2});
  ls.push_back({10, 3});
  ls.push_back({10, 4});
  ls.push_back({10, 5});
  ls.push_back({10, 6});
  ls.push_back({11, 6});
  ls.push_back({12, 6});
  ls.push_back({13, 6});
  ls.push_back({14, 6});
  ls.push_back({15, 6});
  ls.push_back({16, 6});
  ls.push_back({17, 6});

  auto data = full_sim2(os, ls, 0.04, 0.1);
  vector<Scan> scans = data_to_scans(data, ls);
  for (int i = 0; i < scans.size(); i++)
    scans[i].write_scan(i);
}*/

/*void
  look(int argc, char* argv[])
  {
    int nscans = 1;
    if (argc < 7) // nscans,z,N,stride, dist, an tol, data folder
      exit(1);
    nscans = atoi(argv[1]);
    N = atoi(argv[3]);
    angle_tolerance = atof(argv[6]);
    dist_tolerance = atof(argv[5]);
    vector<Scan> scans(nscans);
    for (int i = 0; i < nscans; i++)
      // if (i == nscans - 1) // to look at single scan
      scans[i].read_scan(argv[7], i, atoi(argv[2]));
    cout << "   segmenting first scan" << endl;

    vector<Model*> models;
    Handler h(models);
    // set the default xy to add for odometry
    Vector2d odom_add(-scans[0].odomloc);
    int n = 0;
    for (auto& scan : scans) {
      // scan.write_scan(n++);
      cout << "*********************************************************" <<
  endl; cout << "************************NEW SCAN*************************" <<
  endl; cout << "*********************************************************" <<
  endl;
      // v.add_points(scan.data(all, {0, 1}), "r."); // viz
      Matrix<double, 1, 2> m;
      m << scan.loc(0), scan.loc(1);
      v.add_points(m, "b.");
      m = {
        scans[0].loc(0) + odom_add(1) + scan.odomloc(1),
        scans[0].loc(1) + odom_add(0) + scan.odomloc(0),
      };
      v.add_points(m, "g.");
      vector<Data> preprocessed_scan = h.preprocess_scan(scan);
      int k = 0;
      for (int i = 0; i < preprocessed_scan.size(); i += atoi(argv[4])) {
        printf("data %i:\n", k++);
        VectorXd measurement = preprocessed_scan[i].d;
        cout << "*********************************************************"
             << endl;
        cout << "PROCESSING: " << i << " - " << measurement.transpose() << endl;
        cout << "*********************************************************"
             << endl;
        // if (preprocessed_scan[i].d(2) >= atof(argv[3]) &&
        //     preprocessed_scan[i].d(2) <= atof(argv[4])) // dist req
        v.add_points(Matrix<double, 1, 2>(preprocessed_scan[i].d({ 0, 1 })),
                     "r.");
      }
      cout << "*********************************************************" <<
  endl;
    }
    v.show();
  }*/

void run(int argc, char *argv[]) {
    int nscans = 1;
    if (argc < 8) // nscans,z,N,stride, dist, an tol, env mult, data folder
        exit(1);
    nscans = stoi(argv[1]);
    N = stoi(argv[3]);
    float mult = strtof(argv[7],nullptr);
    angle_tolerance = strtof(argv[6],nullptr);
    dist_tolerance = strtof(argv[5],nullptr);
    vector<Scan> scans(nscans);
    for (int i = 0; i < nscans; i++)
        scans[i].read_scan(string(argv[8]), i, stoi(argv[2]), mult);

    auto *ellipse_model1 = new Model(
            fsn_0, dfsn_0, dop_0, ls_0, ap_ls_0, safety_0, dst_0, init_post_0,
            associate_0, augment_post_0, string(argv[8]) + "/model_e_real.txt");
    /*Model* ellipse_model2 = new Model(fsn_0,
                                      dfsn_0,
                                      dop_0,
                                      ls_2,
                                      ap_ls_2,
                                      safety_2,
                                      dst_0,
                                      init_post_0,
                                      associate_0,
                                      augment_post_0,
                                      string(argv[8]) + "/model_e_real.txt");
    Model* ellipse_model3 = new Model(fsn_0,
                                      dfsn_0,
                                      dop_0,
                                      ls_3,
                                      ap_ls_3,
                                      safety_3,
                                      dst_0,
                                      init_post_0,
                                      associate_0,
                                      augment_post_0,
                                      string(argv[8]) + "/model_e_real.txt");
    Model* line_model = new Model(fsn_1,
                                  dfsn_1,
                                  dop_1,
                                  ls_1,
                                  safety_1,
                                  fsn_1, // same dist func
                                  init_post_1,
                                  associate_1,
                                  augment_post_1,
                                  string(argv[8]) + "/model1_real.txt");*/

    cout << "   segmenting first scan" << endl;
    // e1 0.3 0.7
    // e2 0.8 1.2
    // e3 1.3 1.7
    vector<Model *> models;
    models.push_back(ellipse_model1);
    // models.push_back(ellipse_model2);
    // models.push_back(ellipse_model3);
    //  models.push_back(line_model);
    Handler h(models);
    // v.add_points(scans.back().data(all, {0, 1}), "r."); // viz
    int n = 0;
    int cnt = 0;
    for (auto &scan: scans) {
        Visualizer v(cnt++);
        // scan.write_scan(n++);
        cout << "*********************************************************" << endl;
        cout << "************************NEW SCAN*************************"
             << cnt - 1 << endl;
        cout << "*********************************************************" << endl;
        // v.add_points(scan.data(all, {0, 1}), "r."); // viz
        Matrix<double, 1, 2> m;
        m << scan.loc(0), scan.loc(1);

        // for (int l = 0; l < scan.data.rows(); l++)
        //   v.add_points(Matrix<double, 1, 2>(scan.data(l, { 0, 1 })), "red"); //
        //   viz
        vector<Data> preprocessed_scan = h.preprocess_scan(scan, v);
        int k = 0;
        ofstream odbgplt("debug.gpt");
        odbgplt << "set style data lines" << endl
                << "set term pngcairo size 2000,1000 enhanced font \"Times,12\" "
                << endl
                << "datafile = \"debugdata.dat\"" << endl
                << "set logscale y" << endl
                << R"(plot datafile using 1:2 with lines title 'x', \)" << endl
                << R"(datafile using 1:3 with lines title 'y', \)" << endl
                << R"(datafile using 1:4 with lines title 't', \)" << endl
                << R"(datafile using 1:5 with lines title 'a', \)" << endl
                << R"(datafile using 1:6 with lines title 'b', \)" << endl
                << "datafile using 1:7 with lines title \'e\'" << endl
                << "set key outside" << endl;;

        for (int i = 0; i < preprocessed_scan.size(); i += stoi(argv[4])) {
            printf("data %i:\n", k++);
            VectorXd measurement = preprocessed_scan[i].d;
            cout << "*********************************************************"
                 << endl;
            cout << "PROCESSING: " << i << " - " << measurement.transpose() << endl;
            cout << "*********************************************************"
                 << endl;
            // v.add_points(Matrix<double, 1, 2>(preprocessed_scan[i].d({0, 1})),
            //              "r."); // viz
            h.process_measurement(preprocessed_scan[i], v);
            // debugging
            // plot uncertainty after each data point
        }
        // scan over
        h.end_scan();
        for (auto &e: h.map.entities) {
            if (e.p.size() == 6)
                v.add_ellipse(e.p);
            if (e.p.size() == 2)
                v.add_segment(e.p, e.t);
        }
        v.add_points(m, "blue");
        v.save();
        cout << "*********************************************************" << endl;
        for (auto &e: h.map.entities)
            cout << e.p.transpose() << endl;
        cout << "**************************DONE***************************" << endl;
    }
    /*for (auto &e : h.map.entities) {
      if (e.p.size() == 6)
        v.add_ellipse(e.p, "g-");
      if (e.p.size() == 2)
        v.add_segment(e.p, e.t, "g-");
    }*/
    // v.show();
}

int main(int argc, char *argv[]) {
    run(argc, argv);
    // render plots
    system("parallel -j 24 gnuplot {} \">\" {.}.png ::: *.gpt");
    return 0;
}
