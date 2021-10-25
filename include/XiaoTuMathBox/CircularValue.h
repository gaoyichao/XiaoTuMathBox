#ifndef XTMB_CIRCULARVALUE_H
#define XTMB_CIRCULARVALUE_H

#include <cmath>

namespace xiaotu {
namespace math {

#define CircValRangeDef(_Name, _L, _H, _Z)            \
    const double _Name::lower_bound  = (_L);          \
    const double _Name::upper_bound  = (_H);          \
    const double _Name::zero = (_Z);                  \
    const double _Name::range = ((_H)-(_L));          \
    const double _Name::half_range= ((_H)-(_L)) * 0.5;

#define CircValRangeSpc(_Name, _L, _H, _Z)            \
    struct _Name                                      \
    {                                                 \
        /* [lower_bound, upper_bound) */              \
        static const double lower_bound;              \
        static const double upper_bound;              \
        static const double zero;                     \
        static const double range;                    \
        static const double half_range;               \
                                                      \
        static_assert((_H>_L) && (_Z>=_L) && (_Z<_H), \
                    ": Range not valid");             \
    };

    template <typename RangeType>
    class CircularValue {
        public:
            CircularValue()
                : value(RangeType::zero)
            { }

            CircularValue(double r)
            {
                value = Wrap(r);
            }

            CircularValue(CircularValue const & c)
            {
                value = c.value;
            }

            CircularValue & operator= (double r)
            {
                value = Wrap(r);
                return *this;
            }

            operator double() const
            {
                return value;
            }

        public:
            static double GetLowerBound() { return RangeType::lower_bound; }
            static double GetUpperBound() { return RangeType::upper_bound; }
            static double GetZero() { return RangeType::zero; }
            static double GetRange() { return RangeType::range; }
            static double GetHalfRange() { return RangeType::half_range; }

            static double Wrap(double f)
            {
                while (f > RangeType::upper_bound)
                    f -= RangeType::range;
                while (f < RangeType::lower_bound)
                    f += RangeType::range;
                return f;
            }

            static bool IsInRange(double r)
            {
                return (r >= RangeType::lower_bound && r < RangeType::upper_bound);
            }

            /*
             * ShortestWalk - 从 c1 走到 c2 的最短距离
             */
            static double ShortestWalk(CircularValue const & c1, CircularValue const & c2)
            {
                double d = c2.value - c1.value;
                if (d < -RangeType::half_range)
                    return d + RangeType::range;
                if (d >= RangeType::half_range)
                    return d - RangeType::range;
                return d;
            }

            friend double operator - (CircularValue<RangeType> const & c1, CircularValue<RangeType> const & c2)
            {
                return CircularValue<RangeType>::ShortestWalk(c2, c1);
            }

        private:
            double value;
    };
   
    CircValRangeSpc(SignedRadRange, -M_PI, M_PI, 0.0);

}
}

#endif

