/************************************************************************
 * 
 * 基于拉格朗日多项式的数值积分
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=Lagrange多项式数值积分
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_INTEGRAL_H
#define XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_INTEGRAL_H


namespace xiaotu {

    template <typename Scalar>
    class LagrangePolynomialIntegral
    {
        public:
            LagrangePolynomialIntegral(Scalar a, Scalar b, std::vector<Scalar> const & x)
                : mLowerBound(a), mHigherBound(b)
            {
                assert(x.size() > 1);

                mCs.resize(x.size());
                for (size_t i = 0; i < x.size(); i++) {
                    auto p = LagrangeBasis(x, i).Integral();
                    mCs[i] = p(b) - p(a);
                }
            }

            Scalar Integral(std::vector<Scalar> const & y)
            {
                assert(y.size() == mCs.size());
                Scalar re = 0;
                for (size_t i = 0; i < y.size(); i++) {
                    re += mCs[i] * y[i];
                }
                return re;
            }

        private:
            Scalar mLowerBound;
            Scalar mHigherBound;

            //! 积分系数
            std::vector<Scalar> mCs;
    };
}

#endif