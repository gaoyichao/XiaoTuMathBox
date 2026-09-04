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

/////////////////////////////////////////////////////////////////////////
//
// 一些关于多项式的有用工具
//
/////////////////////////////////////////////////////////////////////////

namespace xiaotu {

    /**
     * @brief 拉格朗日多项式插值
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] y_nodes 采样点 y 坐标
     * @param [in] x 插值点
     */
    template <typename DataType>
    DataType LagrangeInterpolation(std::vector<DataType> const & x_nodes, 
                                   std::vector<DataType> const & y_nodes, 
                                   DataType x)
    {
        size_t n = x_nodes.size();
        assert(y_nodes.size() == n);

        DataType result = 0;
        for (int i = 0; i < n; ++i) {
            DataType Li = 1.0;
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;
                Li *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
            }
            result += y_nodes[i] * Li;
        }

        return result;
    }


    /**
     * @brief 计算拉格朗日插值多项式的系数
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] y_nodes 采样点 y 坐标
     * @return 拉格朗日多项式
     */
    template <typename Polynomial>
    Polynomial LagrangePolynomial(std::vector<typename Polynomial::Scalar> const & x_nodes, 
                                    std::vector<typename Polynomial::Scalar> const & y_nodes)
    {
        using DataType = typename Polynomial::Scalar;

        size_t n = x_nodes.size();
        assert(y_nodes.size() == n);

        Polynomial re;
        re.ReAlloc(n-1);

        for (int i = 0; i < n; ++i) {
            std::vector<DataType> L_i = {1.0}; 
            
            // L_i(x) = \prod_{j \neq i} (x - x_j) / (x_i - x_j)
            DataType denominator = 1.0;
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;

                // L_i = L_i * (x - x_nodes[j])
                std::vector<DataType> next_L_i(L_i.size() + 1, 0.0);
                for (int k = 0; k < L_i.size(); ++k) {
                    next_L_i[k + 1] += L_i[k];
                    next_L_i[k]     -= L_i[k] * x_nodes[j];
                }
                L_i = next_L_i;
                
                denominator *= (x_nodes[i] - x_nodes[j]);
            }
            
            DataType scale = y_nodes[i] / denominator;
            for (int k = 0; k < L_i.size(); ++k) {
                re[k] += scale * L_i[k];
            }
        }

        return re;
    }




}



#endif

