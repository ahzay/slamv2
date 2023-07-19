//
// Created by user on 7/12/23.
//

#include "entitymap.h"

bool EntityMap::associate_data(Data &d) {
    //cout << "associating :" << d.get_xy().transpose() << endl;
    double min_mahalanobis = std::numeric_limits<double>::max();
    Entity *min_ptr = nullptr;
    for (auto &e: entities) {
        if (e.m->_closed_fs) { // check if model is closed
            if (e.m->fs(e, d) <= 0) {
                //cout << "point is inside/on closed curve! ";
                min_mahalanobis = 0;
                min_ptr = &e;
                break;
            }
        }
        // check
        double mahalanobis = e.mahalanobis(d);
        if (e.m->associate(e, d)) {
            if (abs(mahalanobis) <= min_mahalanobis) {     // post req
                min_mahalanobis = abs(mahalanobis);
                min_ptr = &e;
            }
        } else {
            //cout << "post req failed!! " << mahalanobis << endl;
        }
    }
    if (min_ptr) {
        if ((min_mahalanobis <= min_ptr->m->_mahalanobis_aug &&
             min_mahalanobis >= min_ptr->m->_mahalanobis_aug_min) ||
            min_mahalanobis == 0) {
            // cout << "*********ASSOCIATED !!!!!! ************" << min_mahalanobis
            //     << ", " << min_ptr->_p.transpose() << endl;
            d._e = min_ptr;
            d._e->m->augment_post(*d._e, d);
            return true;
        } else {
            // cout << "*********not associated :( ************" << min_mahalanobis
            //     << endl;
            return false;
        }
    }
    return false;
}

void EntityMap::add_augment(const Data &data) {
    //if (data._e == nullptr) // safety (can be removed)
    //    cout << "tried to augment with non-associated data !!!!" << endl;
    for (auto &iekf: iekfs) iekf.add(data);
}

void EntityMap::add_entity(const Entity &e, const Aggregate &a) {
    entities.push_back(e);
    if (!is_psd(entities.back()._E)) { // TODO: maybe cleaner propagation?
        cout << "Entity error matrix is non semi-positive definite, eliminating entity" << endl;
        entities.pop_back();
    } else {
        iekfs.emplace_back(&entities.back());
        entities.back().m->init_post(entities.back(), a);
    }
}

void EntityMap::run_augment() {
    for (auto &iekf: iekfs)
        if (!iekf.a._data_vector.empty())
            iekf.update();
}

EntityMap::EntityMap() = default;




