#pragma once

#include <cmath>
#include <utility>

#include <armadillo>

namespace KF::Math {

// Based on https://stackoverflow.com/a/64247207
template <class S>
constexpr std::pair<S, S> sincos(S arg) {
    return {std::sin(arg), std::cos(arg)};
}

}  // namespace KF::Math
