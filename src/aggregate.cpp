//
// Created by user on 7/10/23.
//

#include "aggregate.h"
#include <boost/graph/simple_point.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>

using namespace boost;

MatrixX2d Aggregate::get_measurement_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i]._measurement;
    return mat;
}

MatrixX2d Aggregate::get_rotated_measurement_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i].get_rotated_measurement();
    return mat;
}

void Aggregate::push_back(const Data &data) {
    _data_vector.push_back(data);
    _data_vector.back().normalize();
    if (!_pose.isApprox(data._pose))
        change_referential(data._pose);
}


void Aggregate::flush(const Data &data) {
    flush();
    push_back(data);

}


void Aggregate::flush() {
    //_data_vector.clear();
    // new flush works with age
    for (auto &d: _data_vector)
        d.life--;
    // COMMENTED FOR TESTING
    erase_if(_data_vector, [](auto d) { return d.life < 1; });
    //_pose.setConstant(NAN);
}

/*Aggregate::Aggregate(const MatrixX2d &mat, const Vector3d &pose) {
    _pose = pose;
    for (int i = 0; i < mat.rows(); i++) {
        _data_vector.emplace_back(mat(i, all), _pose);
    }
}*/

MatrixX2d Aggregate::get_xy_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i].get_xy();
    return mat;
}

void Aggregate::self_sort() {
    sort(_data_vector.begin(), _data_vector.end(),
         [](const auto &d1, const auto &d2) -> bool {
             return d1._measurement(1) < d2._measurement(1);
         });
}

void Aggregate::change_referential(Vector3d new_pose) {
    for (auto &d: _data_vector)
        d.change_referential(new_pose);
    _pose = new_pose;
    // sort
    self_sort();
}

void Aggregate::clear() {
    _data_vector.clear();
    _pose.setConstant(NAN);
}

void Aggregate::push_back(const Aggregate &a) {
    for (const auto &d: a._data_vector)
        push_back(d);
}

void Aggregate::reorder() {
    if (_data_vector.empty()) return;
    self_sort();
    auto max_difference_iter = _data_vector.begin();
    double max_difference = (_data_vector.front().get_xy() - _data_vector.back().get_xy()).norm();

    for (auto it = _data_vector.begin() + 1; it != _data_vector.end(); ++it) {
        double difference = (it->get_xy() - (it - 1)->get_xy()).norm();
        if (difference > max_difference) {
            max_difference = difference;
            max_difference_iter = it;
        }
    }
    rotate(_data_vector.begin(), max_difference_iter, _data_vector.end());
}

typedef adjacency_matrix<undirectedS, no_property, property<edge_weight_t, double> > Graph;
typedef graph_traits<Graph>::vertex_descriptor VertexDescriptor;

double euclideanDistance(const Vector2d &p1, const Vector2d &p2) {
    return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2));
}


/*void Aggregate::self_sort() {
    const int numNodes = _data_vector.size();

    Graph g(numNodes);
    property_map<Graph, edge_weight_t>::type weightMap = get(edge_weight, g);

    // Populate graph with edge weights representing euclidean distances
    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            double distance = euclideanDistance(_data_vector[i].get_xy(), _data_vector[j].get_xy());
            add_edge(i, j, distance, g);
        }
    }

    // Compute approximate TSP tour
    vector<VertexDescriptor> tspTour;
    metric_tsp_approx_tour(g, back_inserter(tspTour));

    // Reorder dataPoints according to TSP tour
    vector<Data> orderedData;
    for (auto v: tspTour) {
        orderedData.push_back(_data_vector[v]);
    }
    _data_vector = orderedData;
}*/


Aggregate::Aggregate() = default;
