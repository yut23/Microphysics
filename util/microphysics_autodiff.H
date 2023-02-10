#ifndef _microphysics_autodiff_H
#define _microphysics_autodiff_H

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <autodiff/forward/dual.hpp>

namespace autodiff {
// add a couple of missing functions
namespace detail {

// natural logarithm of 1+x
struct Log1pOp {};

template<typename R>
using Log1pExpr = UnaryExpr<Log1pOp, R>;

template<typename R, Requires<isExpr<R>> = true> constexpr auto log1p(R&& r) -> Log1pExpr<R> { return { r }; }

template<typename T, typename G>
constexpr void apply(Dual<T, G>& self, Log1pOp)
{
    const T aux = One<T>() / (1.0 + self.val);
    self.val = std::log1p(self.val);
    self.grad *= aux;
}


// cube root
struct CbrtOp {};

template<typename R>
using CbrtExpr = UnaryExpr<CbrtOp, R>;

template <typename R, Requires<isExpr<R>> = true> constexpr auto cbrt(R&& r) -> CbrtExpr<R> { return { r }; }

template<typename T, typename G>
constexpr void apply(Dual<T, G>& self, CbrtOp)
{
    self.val = std::cbrt(self.val);
    self.grad *= 1.0 / (3.0 * self.val * self.val);
}

} // namespace detail

// Redefine dual to use amrex::Real instead of double
using dual = HigherOrderDual<1, amrex::Real>;

// A new namespace that has both the STL math functions and the overloads for
// dual numbers, so we can call the same function whether we're using
// autodiff::dual or amrex::Real.
namespace math_functions {

using std::abs, autodiff::detail::abs;
using std::acos, autodiff::detail::acos;
using std::asin, autodiff::detail::asin;
using std::atan, autodiff::detail::atan;
using std::atan2, autodiff::detail::atan2;
using std::cos, autodiff::detail::cos;
using std::exp, autodiff::detail::exp;
using std::log10, autodiff::detail::log10;
using std::log, autodiff::detail::log;
using std::pow, autodiff::detail::pow;
using std::sin, autodiff::detail::sin;
using std::sqrt, autodiff::detail::sqrt;
using std::tan, autodiff::detail::tan;
using std::cosh, autodiff::detail::cosh;
using std::sinh, autodiff::detail::sinh;
using std::tanh, autodiff::detail::tanh;
using std::erf, autodiff::detail::erf;
using std::hypot, autodiff::detail::hypot;

using std::log1p, autodiff::detail::log1p;
using std::cbrt, autodiff::detail::cbrt;

using std::min, autodiff::detail::min;
using std::max, autodiff::detail::max;

} // namespace math_functions

} // namespace autodiff

#endif
