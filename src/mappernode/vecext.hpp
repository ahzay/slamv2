#ifndef VECEXT_HPP
#define VECEXT_HPP

#include <cmath>
#include <vector>

using namespace std;

// TODO: parallelize loops on vectors
template<typename T>
class vecext : public std::vector<T> {
public:
    using vector<T>::vector; // use the constructors from vector
    // T &operator[](int i) { return vector<T>::at(i); } // range-checked
    // const T &operator[](int i) const { return vector<T>::at(i); } //
    // range-checked
    // minima function
    int minimas_idx(int which) {
        for (int i = 1; i < this->size() - 1; i++) {
            T prev = abs(this->at(i - 1));
            T curr = abs(this->at(i));
            T next = abs(this->at(i + 1));
            if (prev > curr && curr < next)
                return i;
        }
        return -1;
    }

    // sum function
    T sum() {
        T s = 0;
        for (T &d: *this)
            s += d;
        return s;
    }

    // mean function
    T mean() {
        T m = this->sum() / this->size();
        return m;
    }

    // first operations with scalars
    vecext operator+(const T s) const {
        vecext<T> v = *this;
        for (T &d: v)
            d += s;
        return v;
    }

    vecext operator-(const T s) const { return *this + (-s); }

    vecext operator*(const T s) const {
        vecext<T> v = *this;
        for (T &d: v)
            d *= s;
        return v;
    }

    vecext operator/(const T s) const {
        vecext<T> v = *this;
        for (T &d: v)
            d /= s;
        return v;
    }

    //  ops with vectors
    // TODO: ADD SAME-SIZE CHECK
    vecext operator+(const vecext<T> v) const {
        vecext<T> k(*this);
        for (int i = 0; i < k.size(); i++) {
            k[i] = this->at(i) + v.at(i);
        }
        return k;
    }

    vecext operator-(const vecext<T> v) const {
        vecext<T> k(*this);
        for (int i = 0; i < k.size(); i++) {
            k[i] = this->at(i) - v.at(i);
        }
        return k;
    }

    vecext operator*(const vecext<T> v) const {
        vecext<T> k(*this);
        for (int i = 0; i < k.size(); i++) {
            k[i] = this->at(i) * v.at(i);
        }
        return k;
    }
};

template<typename T>
vecext<T> pow(const vecext<T> v, const T s) {
    vecext<T> d(v);
    for (T &j: d)
        j = pow(j, s);
    return d;
}

// inverse operators outside the class
template<typename T>
vecext<T> operator+(T s, vecext<T> v) {
    for (T &d: v)
        d += s;
    return v;
}

template<typename T>
vecext<T> operator-(T s, vecext<T> v) {
    for (T &d: v)
        d -= s;
    return v;
}

template<typename T>
vecext<T> operator*(T s, vecext<T> v) {
    for (T &d: v)
        d *= s;
    return v;
}

template<typename T>
vecext<T> operator/(T s, vecext<T> v) {
    for (T &d: v)
        d /= s;
    return v;
}
// vecmult

// linspace
template<typename T>
vecext<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N - 1);
    vecext<T> xs(N);
    typename vecext<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

// covariance
template<typename T>
T cov(vecext<T> x, vecext<T> y) {
    vecext<T> xm = x - x.mean();
    vecext<T> ym = y - y.mean();
    vecext<T> m = xm * ym;
    T cov = m.sum() / (m.size() - 1);
    return cov;
}

#endif // VECEXT_HPP
