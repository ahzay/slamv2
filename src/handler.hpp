#ifndef HANDLER_HPP
#define HANDLER_HPP
#include "fsm.hpp"

class Entity_map {
public:
  Entity_map() { odbg.open("debugdata.dat"); }
  bool associate_data(Data &data) {
    cout << "associating :" << data.d.transpose() << endl;
    double min_mahalanobis = std::numeric_limits<double>::max();
    Entity *min_ptr = nullptr;
    for (auto &e : entities) {
      if (e.m->_closed_fs) { // check if model is closed
        if (e.m->_fs(e.p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) <= 0) {
          cout << "point is inside/on closed curve! ";
          min_mahalanobis = 0;
          min_ptr = &e;
          break;
        }
      }
      // check
      double mahalanobis = e.mahalanobis(data);
      if (abs(mahalanobis) <= min_mahalanobis && // maha req
          e.m->_associate(e.p, e.t, data)) {     // post req
        min_mahalanobis = abs(mahalanobis);
        min_ptr = &e;
      }
    }
    if (min_ptr) {
      if (min_mahalanobis <= min_ptr->m->_mahalanobis_aug) {
        cout << "*********ASSOCIATED !!!!!! ************" << min_mahalanobis
             << ", " << min_ptr->p.transpose() << endl;
        data.e = min_ptr;
        return 1;
      } else {
        cout << "*********not associated :( ************" << min_mahalanobis
             << endl;
        return 0;
      }
    }
    return 0;
  }

  void add_augment(Data &data) {
    if (data.e == nullptr) // safety (can be removed)
      cout << "tried to augment with non-associated data !!!!" << endl;
    else {
      // Iekf iekf(data.e);
      // iekf.update(data.d, data.p, 0.1);
      // augment with correct iekf
      // probably requires complete restructure to be made efficient
      for (auto &iekf : iekfs)
        if (iekf._e == data.e)
          iekf.add(data.d, data.p);
      // data.e->m->_augment_post(data.e->p, data.e->t, data); // post
    }
  }
  void run_augment() {
    for (auto &iekf : iekfs)
      if (iekf._m0.cols()) {
        iekf.update(0.001);
        odbg << cnt++ << " " << iekf._e->E.diagonal().transpose() << endl;
      }
  }
  // attributes
  vector<Entity> entities; // to be processed map
  vector<Iekf> iekfs;      // init when new entity, reset each scan
  ofstream odbg;
  int cnt = 0;
};

class Handler {
public:
  // non state functions
  Handler(const vector<Model *> models) : _models(models) {}
  vector<Data> preprocess_scan(Scan scan, Visualizer &v) {
    vector<Data> d_v1, d_v2;
    cout << "associating one way!!! rows: " << scan.data.rows() << endl;
    for (int i = 0; i < scan.data.rows(); i++) {
      cout << i << endl;
      VectorXd vbuf = scan.data(i, all);
      Data d(vbuf, scan.loc);
      if (map.associate_data(d)) {
        map.add_augment(d);
        v.add_points(Matrix<double, 1, 2>(d.d({0, 1})), "olive");
      } else {
        v.add_points(Matrix<double, 1, 2>(d.d({0, 1})), "red");
        d_v1.push_back(d);
      }
    }
    // grid for association viz
    /*for (double i = -11; i < 12; i += 0.2)
      for (double j = -20; j < 21; j += 0.2) {
        // transform into d an
        double relx = i - scan.loc(0);
        double rely = j - scan.loc(1);
        double an = atan2(rely, relx);
        double d = sqrt(relx * relx + rely * rely);
        Data dat(Vector<double, 4>({ i, j, d, an }), scan.loc);
        if (map.associate_data(dat))
          v.add_points(Matrix<double, 1, 2>(dat.d({ 0, 1 })), "purple");
        else
          v.add_points(Matrix<double, 1, 2>(dat.d({ 0, 1 })), "yellow");
      }*/

    /*cout << "other way!!!" << endl;
    for (int i = d_v1.size() - 1; i >= 0; i--) {
      if (map.associate_data(d_v1[i]))
        map.augment_from_data(d_v1[i]);
      else
        d_v2.push_back(d_v1[i]); // only add non associated (2 passes)
    }*/
    // return d_v2;
    map.run_augment();
    return d_v1;
  }
  // state functions
  void f_begin() {
    a.flush(_data);
    n = 1;
    state = s_continuous;
    fsms.clear();
  }
  void f_continuous() {
    if (a.check_continuity(_data)) {
      cout << "not continuous!" << endl;
      state = s_begin;
      return;
    }
    a.push_back(_data);
    if (n < N - 1)
      n++;
    else {
      state = s_fsm;
      f_fsm();
    }
  }
  void f_fsm() {
    // if they are empty, initialize them for all models
    if (fsms.empty())
      for (auto m : _models)
        fsms.push_back(Fsm(m, a));
    // instead check if existing to keep order
    // fsms_are_done.clear();
    for (auto &fsm : fsms) {

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
  void f_close_fsms() {
    Fsm *min_fsm = fsms_are_done.back();
    if (min_fsm != nullptr) {
      if (min_fsm->state > 3) { // TODO: make sure we are not adding sink fsms
        // min_fsm->f_least_squares();
        min_fsm->_ekf->_p =
            min_fsm->_model->_ls(min_fsm->_a.pose, min_fsm->_a.get_mat());
        min_fsm->_ekf->_E = min_fsm->_model->_dop(
            min_fsm->_ekf->_p, min_fsm->_a.pose, min_fsm->_a.get_mat(),
            min_fsm->_model->_dop_sigma);
        map.entities.push_back(
            Entity(min_fsm->_a.m, min_fsm->_ekf->_p, min_fsm->_ekf->_E));
        map.entities.back().m->_init_post(map.entities.back().p, // post
                                          map.entities.back().t, min_fsm->_a);
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
  void end_scan() {
    // stl find "minimal" fsm (see < in fsm)
    // first make a vector of only running fsms
    vector<Fsm> running_fsms;
    std::copy_if(
        fsms.begin(), fsms.end(), std::back_inserter(running_fsms),
        [](const Fsm &fsm) -> bool { return (fsm.state == s_fsm_strict); });
    // then least squares
    cout << "running fsms: " << running_fsms.size() << endl;
    for (auto &fsm : running_fsms) {
      // fsm.f_least_squares();
      fsm._ekf->_p = fsm._model->_ls(fsm._a.pose, fsm._a.get_mat());
      fsm._ekf->_E = fsm._model->_dop(fsm._ekf->_p, fsm._a.pose,
                                      fsm._a.get_mat(), fsm._model->_dop_sigma);
      cout << "E: " << endl << fsm._ekf->_E.diagonal().transpose() << endl;
    }

    // then find the minimum
    vector<Fsm>::iterator min_it =
        min_element(running_fsms.begin(), running_fsms.end());
    if (min_it != running_fsms.end()) {
      map.entities.push_back(
          Entity(min_it->_model, min_it->_ekf->_p, min_it->_ekf->_E));
      map.entities.back().m->_init_post(map.entities.back().p, // post
                                        map.entities.back().t, min_it->_a);
    }
    // remove uncertain ellipses
    // std::erase_if(map.entities, [](const auto& e) {
    //  return e.E.diagonal().array().abs().sum() > 100;
    //});
    // reset iekfs for map
    map.iekfs.clear();
    for (auto &e : map.entities)
      map.iekfs.push_back(Iekf(&e));
    // reset attributes
    fsms_are_done.clear();
    fsms.clear();
    n = 0;
    state = s_begin;
  }

  // main functions
  State process_measurement(Data &data, Visualizer &v) {
    //
    cout << "ellipses so far: " << endl;
    for (auto &e : map.entities)
      if (e.m->_model_index != 1)
        cout << "p: " << e.p.transpose() << endl
             << "E: " << e.E.diagonal().transpose() << endl;
    //
    _data = data;
    // try associating and augmenting (since not associated)
    /*if (map.associate_data(data)) {
      // v.add_points(Matrix<double, 1, 2>(data.d(0, { 0, 1 })), "red");
      map.augment_from_data(data);
      state = s_begin;
      return state;
    } else
      v.add_points(Matrix<double, 1, 2>(data.d({ 0, 1 })), "olive");*/
    //
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
  // attributes
  const vector<Model *> _models; // permanent
  // VectorXd _data, _pose;         // buffers, pose should stay during scan
  Data _data;
  unsigned short n = 0;
  State state = s_begin;
  Entity_map map;
  Aggregate a; // buffer aggregate
  vector<Fsm> fsms;
  vector<Fsm *> fsms_are_done;
  ofstream odbg;
};

#endif // HANDLER_HPP
