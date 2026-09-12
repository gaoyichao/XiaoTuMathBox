/************************************************************************
 * 
 * 牛顿差商多项式(Newton Divided Difference)
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=牛顿差商多项式
 * 
 ***********************************************************************/


#ifndef XTMB_NUMERICAL_NEWTON_DIVIDED_DIFFERENCE_H
#define XTMB_NUMERICAL_NEWTON_DIVIDED_DIFFERENCE_H


#include <vector>
#include <functional>
#include <iomanip>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    template <typename Scalar>
    class NewtonDividedDifference
    {
        public:
             /**
             * @brief 构造函数
             * 
             * @param [in] x 采样点 x 列表
             * @param [in] y 对应 x 列表的采样值
             */
            NewtonDividedDifference(
                    std::vector<Scalar> const & x,
                    std::vector<Scalar> const & y)
                : mXs(x)
            {
                assert(!x.empty());
                assert(x.size() == y.size());
        
                InitTable(y);
            }
        
            /**
             * @brief 计算插值
             */
            Scalar Evaluate(Scalar x) const
            {
                int n = mXs.size();
                Scalar result = mTable[0][n - 1]; 
                
                for (int i = n - 2; i >= 0; --i)
                    result = result * (x - mXs[i]) + mTable[0][i];
                return result;
            }

            /**
             * @brief 计算插值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }

            /**
             * @brief 新增采样点并更新差商表
             * 
             * @param [in] x 新增采样点的 x 坐标
             * @param [in] y 新增采样点的 y 坐标
             */
            void AddPoint(Scalar x, Scalar y)
            {
                int n = mXs.size();
                mXs.push_back(x);

                mTable.resize(n + 1);
                for (int i = 0; i <= n; ++i)
                    mTable[i].resize(n + 1, 0.0);
                
                mTable[n][0] = y;
                for (int k = 1; k <= n; ++k) {
                    mTable[n - k][k] = (mTable[n - k + 1][k - 1] - mTable[n - k][k - 1]) 
                                    / (mXs[n] - mXs[n - k]);
                }
            }
            
        private:

            /**
             * @brief 构造通用差商表
             */
            void InitTable(std::vector<Scalar> const & y)
            {
                int n = mXs.size();
                mTable.assign(n, std::vector<Scalar>(n, 0.0));
                for (int i = 0; i < n; ++i)
                    mTable[i][0] = y[i];
        
                for (int j = 1; j < n; ++j) {
                    for (int i = 0; i < n - j; ++i) {
                        mTable[i][j] = (mTable[i+1][j-1] - mTable[i][j-1])
                                     / (mXs[i+j] - mXs[i]);
                    }
                }
            }
    

            friend std::ostream & operator << (
                    std::ostream& os,
                    NewtonDividedDifference const & ndd)
            {
                std::ios_base::fmtflags original_flags = os.flags();
                std::streamsize original_precision = os.precision();

                os << std::fixed << std::setprecision(6);
                os << "\n========================================" << std::endl;

                int n = ndd.mXs.size();
                for (int i = 0; i < n; ++i) {
                    os << "x[" << i << "] = " << std::setw(10) << ndd.mXs[i] << " | ";
                    
                    for (int j = 0; j < (n-i); ++j)
                        os << std::setw(12) << ndd.mTable[i][j] << " ";
                    os << std::endl;
                }
                os << "========================================" << std::endl;

                os.flags(original_flags);
                os.precision(original_precision);
                return os;
            }

        private:
            //! 采样点 x
            std::vector<Scalar> mXs;
            //! 差商表
            std::vector<std::vector<Scalar>> mTable;
    };

}

#endif
