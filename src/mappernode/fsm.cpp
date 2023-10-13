//
// Created by user on 7/12/23.
//

#include "fsm.h"

Fsm::Fsm(Model *model, const Aggregate &a) : _a(a),
                                             e(model, a),
                                             ekf(e),
                                             _data(Vector2d(),
                                                   Vector3d(),
                                                   model->_options.data_longevity) {}


State Fsm::process_measurement(Data &data) {
    _data = data;
    switch (state) {
        case s_fsm_flexible:
            cout << "fsm flexible" << endl;
            f_flexible();
            break;
        case s_fsm_strict:
            cout << "fsm strict" << endl;
            f_strict();
            break;
        case s_fsm_close:
            cout << "fsm close" << endl;
            break;
        case s_fsm_sink:
            cout << "fsm sink" << endl;
            break;
        default:
            break;
    }
    return state;
}

void Fsm::f_flexible() {
    if (ekf.update(_data, false))
        state = s_fsm_sink; // goto sink
    else {
        _a.push_back(_data);
        state = s_fsm_strict; // goto strict next
    }
}

void Fsm::f_strict() {
    if (ekf.update(_data, true)) {
        cout << "going to s_fsm_least_squares" << endl;
        f_least_squares(false); // LS (it is called from here, state changed inside)
        return;
    }
    // else add data and stay in strict
    _a.push_back(_data);
}

void Fsm::f_least_squares(const bool is_fsm_end) {
    // here we LS, DOP and reset EKF
    e.m->ls(e, _a, 0);
    if (!is_fsm_end)
        if (ekf.update(_data, true)) {
            cout << "FSM CLOSING!!!!" << endl;
            state = s_fsm_close; // goto close next, what to do with data
        } else {
            _a.push_back(_data);
            state = s_fsm_strict; // goto strict next
        }
}

/*bool Fsm::operator<(const Fsm &other_fsm) const {
    // assuming we are in flexible (TODO: account if not?)
    return (ekf._mahalanobis < other_fsm.ekf._mahalanobis);
}*/