//
// Created by user on 7/12/23.
//

#include "handler.h"

Handler::Handler(const deque<Model *> &models,
                 const CmdLineOptions &options) : _models(models),
                                                  _data(Data(Vector2d(),
                                                             Vector3d())),
                                                  o(options), ps(options), ua(options),
                                                  v(nullptr) {}

void Handler::preprocess_scan(Scan &scan) {
    cout << "associating" << endl;
    Aggregate aa; // associated aggregate (only for display)
    for (auto &d: scan._data_vector) {
        if (map.associate_data(d)) {
            map.add_augment(d);
            aa.push_back(d);
            //v.add_points(Matrix<double, 1, 2>(d.d({0, 1})), "olive");
        } else {
            //v.add_points(Matrix<double, 1, 2>(d.d({0, 1})), "red");
            ps.push_back(d);
            //d_v1.push_back(d);
        }
    }
    v->add_aggregate(aa, "green");
    v->add_aggregate(ps, "red");

    // debug grid for association viz
    /*for (double i = -11; i < 12; i += 0.2)
      for (double j = -20; j < 21; j += 0.2) {
        // transform into d an
        double relx = i - scan.loc(0);
        double rely = j - scan.loc(1);
        double an = atan2(rely, relx);
        double d = sqrt(relx * relx + rely * rely);
        Data dat(Vector<double, 4>({i, j, d, an}), scan.loc);
        if (map.associate_data(dat))
          v.add_points(Matrix<double, 1, 2>(dat.d({0, 1})), "purple");
        else
          v.add_points(Matrix<double, 1, 2>(dat.d({0, 1})), "yellow");
      }*/
    map.run_augment();
}

void Handler::process_scan() {
    for (auto &data: ps._data_vector)
        process_measurement(data);
    end_scan();
}

State Handler::process_measurement(Data &data) {
    _data = data;
    switch (state) {
        case s_begin:
            cout << "begin nsm" << endl;
            f_begin();
            break;
        case s_continuous:
            cout << "nsm continuous, n=" << n << endl;
            f_continuous();
            break;
        case s_fsm:
            cout << "nsm fsm, n=" << n << endl;
            f_fsm();
            break;
        default:;
    }
    return state;
}

void Handler::f_begin() {
    // reset
    ua.flush(_data);
    n = 1;
    state = s_continuous;
    fsms.clear();
}

void Handler::f_continuous() {
    if (ua.check_continuity(_data)) {
        cout << "not continuous!" << endl;
        state = s_begin;
        return;
    }
    ua.push_back(_data);
    if (n < o.init_npoints - 1)
        n++;
    else {
        state = s_fsm;
        f_fsm();
    }
}

void Handler::f_fsm() {
    // if they are empty, initialize them for all models
    if (fsms.empty())
        for (auto m: _models)
            fsms.emplace_back(m, ua);
    // instead check if existing to keep order
    // fsms_are_done.clear();
    for (auto &fsm: fsms) {
        State buf_state = fsm.process_measurement(_data);
        if (buf_state == s_fsm_close || buf_state == s_fsm_sink) {
            // check if not in fsms_done
            if (find(fsms_are_done.begin(), fsms_are_done.end(), &fsm) ==
                fsms_are_done.end())
                fsms_are_done.push_back(&fsm); // push back if not
        }
        // viz for debugging!
        /*else {
          if (fsm._model->_model_index == 0)
            v.add_ellipse(fsm._ekf->_p, "b-");
          if (fsm._model->_model_index == 1)
            v.add_line(fsm._ekf->_p, "b-");
        }*/
        //
    }
    if (fsms_are_done.size() == fsms.size()) // not all done
        f_close_fsms(); // we assume that the most accurate model has lowest r
}

void Handler::f_close_fsms() {
    Fsm *min_fsm = fsms_are_done.back();
    if (min_fsm != nullptr) {
        if (min_fsm->state > 3) { // TODO: make sure we are not adding sink fsms
            min_fsm->f_least_squares(true);
            map.add_entity(min_fsm->e);
            map.entities.back().m->init_post(map.entities.back(),
                                             min_fsm->_a);
        }
    }

    // reset attributes
    fsms_are_done.clear();
    fsms.clear();
    n = 0;
    state = s_begin;
    // call f_begin as to not lose data point for new cycle
    f_begin();
}

void Handler::end_scan() {
    // stl find "minimal" fsm (see < in fsm)
    // first make a vector of only running fsms
    cout << "Force ending scan: " << endl;
    vector<Fsm> running_fsms;
    std::copy_if(
            fsms.begin(),
            fsms.end(),
            std::back_inserter(running_fsms),
            [](const Fsm &fsm) -> bool { return (fsm.state == s_fsm_strict); });
    // then least squares
    cout << "running fsms: " << running_fsms.size() << endl;
    for (auto &fsm: running_fsms) {
        fsm.f_least_squares(true);
        //cout << "E: " << endl << fsm._ekf->_E.diagonal().transpose() << endl;
    }

    // then find the minimum
    vector<Fsm>::iterator min_it =
            min_element(running_fsms.begin(), running_fsms.end());
    if (min_it != running_fsms.end()) {
        map.add_entity(min_it->e);
        map.entities.back().m->init_post(map.entities.back(),
                                         min_it->_a);
    }
    // remove uncertain ellipses
    // std::erase_if(map.entities, [](const auto& e) {
    //  return e.E.diagonal().array().abs().sum() > 100;
    //});
    // reset iekfs for map
    /*
    map.iekfs.clear();
    for (auto &e: map.entities)
        map.iekfs.push_back(Iekf(&e));*/
    // reset attributes
    // maybe there is a more elegant way for this?
    for (auto &iekf: map.iekfs)
        iekf.a.flush();
    ua.flush();
    ps.flush();
    fsms_are_done.clear();
    fsms.clear();
    n = 0;
    state = s_begin;
    cout << "Force ended scan with entities: " << endl;
    for (const auto &e: map.entities)
        e.print();
}