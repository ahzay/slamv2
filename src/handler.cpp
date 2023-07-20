//
// Created by user on 7/12/23.
//

#include "handler.h"

Handler::Handler(const deque<Model *> &models,
                 const CmdLineOptions &options) : _models(models),
                                                  _data(Data(Vector2d(),
                                                             Vector3d(),
                                                             options.data_longevity)),
                                                  o(options), ps(options),
                                                  ua(options), uua(options), v(nullptr) {}

void Handler::preprocess_scan(Scan &scan) {
    cout << "start preprocessing" << endl;
    // for speedup
    if (uua._data_vector.size() > o.init_npoints * o.npoints_mult) {
        auto tmp = Scan(o);
        mt19937 gen(random_device{}());
        tmp._data_vector.resize(o.init_npoints * o.npoints_mult,
                                Data(Vector2d(),
                                     Vector3d(),
                                     o.data_longevity));
        sample(uua._data_vector.begin(),
               uua._data_vector.end(),
               tmp._data_vector.begin(),
               o.init_npoints * o.npoints_mult, gen);
        uua = tmp;
    }
    //
    cout << "associating" << endl;
    Aggregate aa; // associated aggregate (only for display)
    Scan buf(o), buff(o);
    buf.push_back(uua);
    buf.push_back(scan);
    int prev_size = 0;
    do {
        prev_size = buff._data_vector.size();
        buff.clear();
        for (auto &d: buf._data_vector) {
            if (map.associate_data(d)) {
                map.add_augment(d);
                aa.push_back(d);
                d.life = 0;
            } else {
                buff.push_back(d);
            }
            buf.clear();
            buf.push_back(buff);
        }
    } while (prev_size != buff._data_vector.size());
    ps.push_back(buff);
    /*for (auto &d: uua._data_vector) {
        if (map.associate_data(d)) {
            map.add_augment(d);
            aa.push_back(d);
            d.life = 0;
        } else {
            ps.push_back(d);
        }
    }

    for (auto &d: scan._data_vector) {
        if (map.associate_data(d)) {
            map.add_augment(d);
            aa.push_back(d);
        } else {
            ps.push_back(d);
        }
    }*/
    uua.clear();
    ps.reorder();
    v->add_aggregate(aa, "green");
    v->add_aggregate(ps, "red");
    v->add_data(ps._data_vector.back(), "orange");
    v->add_data(ps._data_vector.front(), "yellow");
    map.run_augment();
    cout << "done preprocessing!" << endl;
}

void Handler::process_scan() {
    for (auto &data: ps._data_vector)
        process_measurement(data);
    end(true);
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
    //v->add_data(_data, "purple");
    adjusted_init_npoints = o.init_npoints / _data._measurement(0) * 3;
    cout << "adjusted init points: " << adjusted_init_npoints << endl;
    ua.clear();
    ua.push_back(_data);
    n = 1;
    state = s_continuous;
    fsms.clear();
}

void Handler::f_continuous() {
    if (ua.check_continuity(_data)) {
        cout << "not continuous!" << endl;
        uua.push_back(ua);
        state = s_begin;
        f_begin();
        return;
    }
    ua.push_back(_data);
    if (n < adjusted_init_npoints - 1)
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
    }
    if (fsms_are_done.size() == fsms.size()) // all done
        end(false); // we assume that the most accurate model has lowest r
}

void Handler::end(bool at_scan_end) {
    // stl find "minimal" fsm (see < in fsm)
    // first make a vector of only running fsms
    cout << "Ending fsms: " << endl;
    for (auto &fsm: fsms)
        fsm.f_least_squares(true);
    vector<Fsm> running_fsms, non_running_fsms;
    partition_copy(
            fsms.begin(),
            fsms.end(),
            back_inserter(running_fsms),
            back_inserter(non_running_fsms),
            [](const Fsm &fsm) -> bool {
                return fsm.state == s_fsm_strict ||
                       fsm.state == s_fsm_close;
            });
    // then least squares
    cout << "running fsms: " << running_fsms.size() << endl;


    // find fsm that has most points inside
    if (!running_fsms.empty()) {
        vector<Fsm>::const_iterator min_it = running_fsms.end();
        if (!at_scan_end) {
            auto most_points_its = all_max_elements(running_fsms,
                                                    [](const Fsm &f1, const Fsm &f2) {
                                                        return f1._a._data_vector.size() < f2._a._data_vector.size();
                                                    });
            if (most_points_its.size() > 1)
                min_it = *min_element(most_points_its.begin(), most_points_its.end(),
                                      [](const auto &f1, const auto &f2) {
                                          return f1->ekf.mahalanobis(f1->_a) <
                                                 f2->ekf.mahalanobis(f2->_a);
                                      });
            else min_it = most_points_its.back();
        }

        if (at_scan_end)
            min_it = min_element(running_fsms.begin(), running_fsms.end(),
                                 [](const Fsm &f1, const Fsm &f2) {
                                     // TODO: move maha fun to FSM
                                     return f1.ekf.mahalanobis(f1._a) <
                                            f2.ekf.mahalanobis(f2._a);
                                 });
        if (min_it != running_fsms.end())
            map.add_entity(min_it->e, min_it->_a);
    }
    if (!non_running_fsms.empty() && !at_scan_end)
        uua.push_back(non_running_fsms.back()._a);
    // TODO: remove uncertain ellipses maybe?
    fsms_are_done.clear();
    fsms.clear();
    state = s_begin;
    if (at_scan_end) {
        uua.push_back(ua); // important
        for (auto &iekf: map.iekfs)
            iekf.a.clear();
        ua.clear();
        ps.clear();
        uua.flush();
        cout << "Force ended scan with entities: " << endl;
        for (const auto &e: map.entities)
            e.print();
    } else f_begin();
}