#ifndef MYMOMENTS_H
#define MYMOMENTS_H

#include "Matrix.h"
#include "myoperators.h"
#include <string>
#include <vector>

template <typename...> class moments;

enum MOMENTS {
    COUNT,
    SUM,
    SUM_SQR,
    MEAN,
    VARIANCE,
    CHI_VALUE,
    T_VALUE,
    T2_VALUE
};

const std::string MOMENTS_string[] = {"count",   "sum",      "sum_sqr",
                                      "mean",    "variance", "chi_value",
                                      "t_value", "T2_value"};

constexpr MOMENTS MOMENTS_VALUES[] = {COUNT,    SUM,       SUM_SQR, MEAN,
                                      VARIANCE, CHI_VALUE, T_VALUE, T2_VALUE};

template <> class moments<double> {

public:
    typedef moments self_type;
    constexpr static auto className =
        my_static_string("moments_") + my_trait<double>::className;
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "count", &self_type::count),
            grammar::field(C<self_type>{}, "sum", &self_type::sum),
            grammar::field(C<self_type>{}, "sum_sqr", &self_type::sum_sqr));
    }

    moments(std::size_t count, double sum, double sum_sqr)
        : count_(count), sum_{sum}, sum_sqr_{sum_sqr} {};

    moments(const std::vector<double> &data)
        : count_(data.size()), sum_(op::sum(data)),
          sum_sqr_(op::sum_sqr(data)) {}

    template <class C, class F>
    moments(const std::vector<C> &data, const F &f)
        : count_(data.size()), sum_(op::sum(data, f)),
          sum_sqr_(op::sum_sqr(data, f)) {}

    std::size_t data_size() const { return 1; }

    std::size_t data_dim(MOMENTS) const { return 0; }

    std::tuple<double> data_row(MOMENTS m) const {
        switch (m) {
        case COUNT:
            return count();
        case SUM:
            return sum();
        case SUM_SQR:
            return sum_sqr();
        case MEAN:
            return mean();
        case VARIANCE:
            return variance();
        case CHI_VALUE:
            return chi_value();
        case T_VALUE:
            return t_value();
        case T2_VALUE:
        default:
            return T2_value();
        }
    }

    std::size_t count() const { return count_; }
    double sum() const { return sum_; }
    double sum_sqr() const { return sum_sqr_; }
    double mean() const { return sum() * (1.0 / count()); }
    double variance() const {
        return sum_sqr() * (1.0 / count()) - sqr(mean());
    }
    double chi_value() const { return mean() / sqrt(variance()); }
    double t_value() const { return mean() / sqrt(variance() / count()); }
    double T2_value() const { return sqr(t_value()); }

    moments() = default;

    moments &operator+=(const moments &other) {
        count_ += other.count();
        sum_ += other.sum();
        sum_sqr_ += other.sum_sqr();
        return *this;
    }

    void push_back(double val) {
        count_++;
        sum_ += val;
        sum_sqr_ += sqr(val);
    }
    moments<double> operator*(double t) {
        double newcount = count() * t;
        double newsum = sum() * t;
        double newsumsqr = sum_sqr() * t;
        return moments<double>(newcount, newsum, newsumsqr);
    }

private:
    double count_;
    double sum_;
    double sum_sqr_;
};

template <typename... T> class moments<std::tuple<T...>> {
    std::tuple<moments<T>...> m_;

public:
    typedef moments self_type;

    constexpr static auto className =
        my_static_string("moments_") + my_trait<std::tuple<T...>>::className;
    static auto get_constructor_fields() {
        return std::tuple_cat(T::get_contructor_fields()...);
    }

    template <std::size_t... Is>
    static std::tuple<moments<T>...>
    build_tuple(std::index_sequence<Is...>,
                const std::vector<std::tuple<T...>> &l) {
        return std::tuple<moments<T>...>(moments<T>{l, &std::get<Is>}...);
    }

    moments(const std::vector<std::tuple<T...>> &l)
        : m_{build_tuple(std::index_sequence_for<T...>{}, l)} {}


    template <class U, class Op, typename... Ts,std::size_t... Is>
    static std::tuple<moments<T>...>
    build_tuple(std::index_sequence<Is...>,
                const std::vector<U> &l, const Op &u, Ts... t) {
        return std::tuple<moments<T>...>(moments<T>{l, [&u, &t...](const U &d) { return get<Is>(u(d, t...)); }}...);
    }



    template <class U, class Op, typename... Ts>
    moments(const std::vector<U> &l, const Op &u, Ts... t)
        : m_{build_tuple(std::index_sequence_for<T...>{}, l,u,t...)}{}

    moments() = default;

    auto data_row(MOMENTS m) const {
        return std::apply(
            [m](auto... x) { return std::tuple_cat(x.data_row(m)...); }, m_);
    }
};

template <typename T> class moments<M_Matrix<T>> {

public:
    typedef moments self_type;
    constexpr static auto className =
        my_static_string("moments_") + my_trait<M_Matrix<T>>::className;
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "count", &self_type::count),
            grammar::field(C<self_type>{}, "sum", &self_type::sum),
            grammar::field(C<self_type>{}, "sum_sqr", &self_type::sum_sqr));
    }

    moments(std::size_t count, M_Matrix<T> sum, M_Matrix<T> sum_sqr)
        : count_(count), sum_{sum}, sum_sqr_{sum_sqr}, mean_{calc_mean()},
          variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {};

    moments(const std::vector<M_Matrix<T>> &data)
        : count_(data.size()), sum_(op::sum(data)), sum_sqr_(op::sum_sqr(data)),
          mean_{calc_mean()}, variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {}

    template <class C, class F>
    moments(const std::vector<C> &data, const F &f)
        : count_(data.size()), sum_(op::sum(data, f)),
          sum_sqr_(op::sum_sqr(data, f)), mean_{calc_mean()},
          variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {}

    std::size_t data_size() const { return size(); }

    std::size_t data_dim(MOMENTS m) const {
        switch (m) {
        case COUNT:
            return 0;
        case SUM:
            return 1;
        case SUM_SQR:
            return 2;
        case MEAN:
            return 1;
        case VARIANCE:
            return 2;
        case CHI_VALUE:
            return 1;
        case T_VALUE:
            return 1;
        case T2_VALUE:
        default:
            return 0;
        }
    }
    std::tuple<double> data_row(MOMENTS m, std::size_t i, std::size_t j) const {
        switch (m) {
        case COUNT:
            return count();
        case SUM:
            return sum()[i];
        case SUM_SQR:
            return sum_sqr()(i, j);
        case MEAN:
            return mean()[i];
        case VARIANCE:
            return variance()(i, j);
        case CHI_VALUE:
            return chi_value()[i];
        case T_VALUE:
            return t_value()[i];
        case T2_VALUE:
        default:
            return T2_value();
        }
    }

    auto count() const { return count_; }
    auto &sum() const { return sum_; }
    auto &sum_sqr() const { return sum_sqr_; }
    auto &mean() const { return mean_; }
    auto &variance() const { return variance_; }

    auto &chi_value() const { return chi_value_; }

    auto &t_value() const { return t_value_; }

    double T2_value() const {
        auto covinv = inv(variance() / count());
        if (!covinv)
            return std::numeric_limits<double>::quiet_NaN();
        else
            return xTSigmaX(mean(), covinv.value());
    }

    moments() = default;

    std::size_t size() const { return sum_.size(); }
    moments &operator+=(const moments &other) {
        count_ += other.count();
        sum_ += other.sum();
        sum_sqr_ += other.sum_sqr();
        mean_ = calc_mean();
        variance_ = calc_variance();
        chi_value_ = calc_chi_value();
        t_value_ = calc_t_value();
    }

private:
    std::size_t count_;
    M_Matrix<T> sum_;
    M_Matrix<T> sum_sqr_;
    M_Matrix<T> mean_;
    M_Matrix<T> variance_;
    M_Matrix<T> chi_value_;
    M_Matrix<T> t_value_;

    M_Matrix<T> calc_variance() {
        return sum_sqr() * (1.0 / count()) - op::sqr(mean());
    }
    M_Matrix<T> calc_mean() { return sum() * (1.0 / count()); }
    M_Matrix<T> calc_chi_value() {
        return mean() * (diag(variance())).apply([](double x) {
            return 1.0 / std::sqrt(x);
        });
    }
    M_Matrix<T> calc_t_value() {
        return mean() * (diag(variance()) / count()).apply([](double x) {
            return 1.0 / std::sqrt(x);
        });
    }
};

template <typename T> class moments_matrix {

public:
    typedef moments_matrix self_type;
    constexpr static auto className =
        my_static_string("moments_matrix_") + my_trait<M_Matrix<T>>::className;
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "count", &self_type::count),
            grammar::field(C<self_type>{}, "sum", &self_type::sum),
            grammar::field(C<self_type>{}, "sum_sqr",
                           &self_type::sum_elem_sqr));
    }

    moments_matrix(std::size_t count, M_Matrix<T> sum, M_Matrix<T> sum_elem_sqr)
        : count_(count), sum_{sum}, sum_elem_sqr_{sum_elem_sqr},
          mean_{calc_mean()}, variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {};

    moments_matrix(const std::vector<M_Matrix<T>> &data)
        : count_(data.size()), sum_(op::sum(data)),
          sum_elem_sqr_(op::sum_sqr_elem(data)), mean_{calc_mean()},
          variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {}

    template <class C, class F>
    moments_matrix(const std::vector<C> &data, const F &f)
        : count_(data.size()), sum_(op::sum(data, f)),
          sum_elem_sqr_(op::sum_sqr_elem(data, f)), mean_{calc_mean()},
          variance_{calc_variance()},
          chi_value_{calc_chi_value()}, t_value_{calc_t_value()} {}

    std::size_t data_size() const { return size(); }

    std::size_t data_dim(MOMENTS m) const {
        switch (m) {
        case COUNT:
            return 0;
        case SUM:
            return 2;
        case SUM_SQR:
            return 2;
        case MEAN:
            return 2;
        case VARIANCE:
            return 2;
        case CHI_VALUE:
            return 2;
        case T_VALUE:
            return 2;
        case T2_VALUE:
        default:
            return 0;
        }
    }
    std::tuple<double> data_row(MOMENTS m, std::size_t i, std::size_t j) const {
        switch (m) {
        case COUNT:
            return count();
        case SUM:
            return sum()(i, j);
        case SUM_SQR:
            return sum_elem_sqr()(i, j);
        case MEAN:
            return mean()(i, j);
        case VARIANCE:
            return variance()(i, j);
        case CHI_VALUE:
            return chi_value()(i, j);
        case T_VALUE:
            return t_value()(i, j);
        case T2_VALUE:
        default:
            return T2_value();
        }
    }

    auto count() const { return count_; }
    auto &sum() const { return sum_; }
    auto &sum_elem_sqr() const { return sum_elem_sqr_; }
    auto &mean() const { return mean_; }
    auto &variance() const { return variance_; }

    auto &chi_value() const { return chi_value_; }

    auto &t_value() const { return t_value_; }

    double T2_value() const { return fullSum(elemMult(t_value(), t_value())); }

    moments_matrix() = default;

    std::size_t size() const { return sum_.size(); }
    moments_matrix &operator+=(const moments_matrix &other) {
        count_ += other.count();
        sum_ += other.sum();
        sum_elem_sqr_ += other.sum_sqr();
        mean_ = calc_mean();
        variance_ = calc_variance();
        chi_value_ = calc_chi_value();
        t_value_ = calc_t_value();
    }

private:
    std::size_t count_;
    M_Matrix<T> sum_;
    M_Matrix<T> sum_elem_sqr_;
    M_Matrix<T> mean_;
    M_Matrix<T> variance_;
    M_Matrix<T> chi_value_;
    M_Matrix<T> t_value_;

    M_Matrix<T> calc_variance() {
        return sum_elem_sqr() * (1.0 / count()) - elemMult(mean(), mean());
    }
    M_Matrix<T> calc_mean() { return sum() * (1.0 / count()); }
    M_Matrix<T> calc_chi_value() {
        return elemDiv(mean(),
                       (variance()).apply([](auto x) { return std::sqrt(x); }));
    }
    M_Matrix<T> calc_t_value() {
        return elemDiv(mean(), (variance() / count()).apply([](auto x) {
            return std::sqrt(x);
        }));
    }
};

template <class T, class = void> struct getMoment { typedef moments<T> type; };

template <class T> struct getMoment<T, std::void_t<typename T::moment_type>> {
    typedef typename T::moment_type type;
};

template <class C> using moments_t = typename getMoment<C>::type;

template <class X, class F, typename... args>
auto moments_by_row(const std::vector<X> &x, const F &f, args... t) {
    typedef std::invoke_result_t<F, X, args...> C;
    typedef typename std::decay_t<C>::value_type T;
    auto &c = std::invoke(f, x[0], t...);

    std::vector<moments_t<T>> out(c.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = moments_t<T>(
            x, [&f, &i, &t...](auto &e) { return std::invoke(f, e, t...)[i]; });
    }
    return out;
}

template <class X, class F, class G>
auto moments_by_row_2(const std::vector<X> &x, const F &f, const G &g) {
    typedef std::invoke_result_t<F, X> C;
    typedef typename std::decay_t<C>::value_type T;
    typedef std::decay_t<std::invoke_result_t<G, T>> K;
    auto &c = std::invoke(f, x[0]);

    std::vector<moments_t<K>> out(c.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = moments_t<K>(x, [&f, &g, &i](auto &e) {
            return std::invoke(g, std::invoke(f, e)[i]);
        });
    }
    return out;
}

#endif // MYMOMENTS_H
