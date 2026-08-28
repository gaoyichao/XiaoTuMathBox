#ifndef XTMB_NUMERICAL_POLYNOMIAL_OPERATORS_H
#define XTMB_NUMERICAL_POLYNOMIAL_OPERATORS_H

#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>

/////////////////////////////////////////////////////////////////////////
//
// 多项式加法
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {


    template <typename Polynomial>
    Polynomial operator + (Polynomial const & a, Polynomial const & b)
    {
        Polynomial re;
        size_t degree = std::max(a.Degree(), b.Degree());
    
        re.ReAlloc(degree);
        for (size_t i = 0; i <= degree; ++i) {
            re[i] = a.Coeff(i) + b.Coeff(i);
        }
        
        re.Normalize();
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator + (typename Polynomial::Scalar const & a, Polynomial const & b)
    {
        Polynomial re = b;
        re[0] += a;
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator + (Polynomial const & a, typename Polynomial::Scalar const & b)
    {
        return b + a;
    }
}

/////////////////////////////////////////////////////////////////////////
//
// 多项式乘法
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    template <typename Polynomial>
    Polynomial operator * (Polynomial const & a, Polynomial const & b)
    {
        if (a.IsZero() || b.IsZero())
            return Polynomial::Zero();
    
        size_t degree = a.Degree() + b.Degree();
        Polynomial re;
    
        re.ReAlloc(degree);
        for (size_t i = 0; i <= a.Degree(); ++i) {
            for (size_t j = 0; j <= b.Degree(); ++j) {
                re[i + j] += a[i] * b[j];
            }
        }
    
        re.Normalize();
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator * (typename Polynomial::Scalar const & a, Polynomial const & b)
    {
        if (Traits<Polynomial>::Abs(a) < Traits<Polynomial>::small_value || b.IsZero())
            return Polynomial::Zero();
    
        Polynomial re = b;
        for (size_t i = 0; i <= b.Degree(); ++i)
            re[i] *= a;
    
        re.Normalize();
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator * (Polynomial const & a, typename Polynomial::Scalar const & b)
    {
        return b * a;
    }    
}

/////////////////////////////////////////////////////////////////////////
//
// 多项式减法
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    template <typename Polynomial>
    Polynomial operator - (Polynomial const & a, Polynomial const & b)
    {
        Polynomial re;
        size_t degree = std::max(a.Degree(), b.Degree());
    
        re.ReAlloc(degree);
        for (size_t i = 0; i <= degree; ++i) {
            re[i] = a.Coeff(i) - b.Coeff(i);
        }
        
        re.Normalize();
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator - (typename Polynomial::Scalar const & a, Polynomial const & b)
    {
        Polynomial re = -1 * b;
        re[0] += a;
        return re;
    }
    
    template <typename Polynomial>
    Polynomial operator - (Polynomial const & a, typename Polynomial::Scalar const & b)
    {
        Polynomial re = b;
        re[0] += a;
        return re;
    }

}

/////////////////////////////////////////////////////////////////////////
//
// 多项式相等判定
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    template <typename Polynomial>
    bool operator == (Polynomial const & a, Polynomial const & b)
    {
        size_t deg = std::max(a.Degree(), b.Degree());
        
        for (size_t i = 0; i <= deg; ++i) {
            if (std::abs(a.Coeff(i) - b.Coeff(i)) >= Traits<Polynomial>::small_value)
                return false;
        }
    
        return true;
    }
    
    template <typename Polynomial>
    bool operator != (Polynomial const & a, Polynomial const & b)
    {
        return !(a == b);
    }
    
    template <typename Polynomial>
    bool operator == (typename Polynomial::Scalar const & a, Polynomial const & b)
    {
        if (b.Degree() > 0)
            return false;
        if (std::abs(a - b[0]) >= Traits<Polynomial>::small_value)
            return false;
        return true;
    }
    
    template <typename Polynomial>
    bool operator != (typename Polynomial::Scalar const & a, Polynomial const & b)
    {
        return !(a == b);
    }
    
    template <typename Polynomial>
    bool operator == (Polynomial const & a, typename Polynomial::Scalar const & b)
    {
        return b == a;
    }
    
    template <typename Polynomial>
    bool operator != (Polynomial const & a, typename Polynomial::Scalar const & b)
    {
        return !(b == a);
    }

}

#endif

