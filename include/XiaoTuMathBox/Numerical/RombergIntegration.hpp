/************************************************************************
 * 
 * Romberg 积分
 * 对梯形法则的 Richardson 外展。
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=Romberg积分
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_ROMBERG_INTEGRATION
#define XTMB_NUMERICAL_ROMBERG_INTEGRATION

#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>


namespace xiaotu {

    template <typename Scalar>
    class RombergIntegration
    {
        public:

            /**
             * @brief 构造函数
             * 
             * @param [in] f 被积函数
             * @param [in] a 积分区间 [a,b]
             * @param [in] b a <= b
             */
            RombergIntegration(std::function<Scalar(Scalar)> && f, Scalar a, Scalar b)
                : _f(f), _a(a), _b(b)
            {
                assert(a < b);

                _k = 0;
                _h = _b - _a;
                mPrevRow.push_back(0.5 * _h * (_f(_a) + _f(_b)));
            }

            /**
             * @brief 获取当前积分值
             */
            Scalar Integral()
            {
                return mPrevRow[_k];
            }

            /**
             * @brief 执行一次 Richardson 外展
             */
            Scalar Extra()
            {
                _k++;
                _h = 0.5 * _h;

                // 只计算当前外展新增采样点
                int num_new_points = 1 << (_k - 1); // 2^(_k - 1)
                Scalar sum_new_f = static_cast<Scalar>(0);
                for (int i = 0; i < num_new_points; ++i) {
                    Scalar new_point = _a + (2 * i + 1) * _h;
                    sum_new_f += _f(new_point);
                }

                mCurrRow.resize(_k+1);
                mCurrRow[0] = 
                // R(k, 0) = _h [f(a)/2 + \sum f(a + i _h) + f(b)/2]
                // R(k, j) = \frac{4^j \cdot R(k,j-1) - R(k-1, j-1)}{4^j-1}
                mCurrRow[0] = 0.5 * mPrevRow[0] + _h * sum_new_f;
                Scalar factor = 4; 
                for (int j = 1; j <= _k; ++j) {
                    mCurrRow[j] = (factor * mCurrRow[j - 1] - mPrevRow[j - 1]) / (factor - 1);
                    factor *= 4; 
                }

                mPrevRow = mCurrRow;
                return mPrevRow[_k];
            }

            /**
             * @brief 自适应 Richardson 外展
             * 
             * @param [in] tol 容忍度, 
             * @param [in] max_iter 最大外展次数
             */
            Scalar Evaluate(Scalar tol = SMALL_VALUE, size_t max_iter = 8)
            {
                Scalar prev = Integral();
                for (size_t i = 0; i < max_iter; ++i) {
                    Scalar curr = Extra();
                    if (std::abs(curr - prev) < tol)
                        return curr;
                }
                return Integral();
            }

            /**
             * @brief 获取当前外展次数
             */
            int NumExtras() const { return _k; }

        private:
            //! 被积函数
            std::function<Scalar(Scalar)> _f;
            //! 积分区间 [_a, _b], _a < _b
            Scalar _a;
            Scalar _b;
            //! 外展次数
            int _k;
            //! 当前积分步长
            Scalar _h;
            // Romberg 表中的上一列和当前列
            std::vector<Scalar> mPrevRow;
            std::vector<Scalar> mCurrRow;
    };


}

#endif
