/************************************************************************
 * 
 * Hermite 插值多项式的差商形式
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=Hermite插值多项式
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_HERMITE_DIVIDED_DIFFERENCE_H
#define XTMB_NUMERICAL_HERMITE_DIVIDED_DIFFERENCE_H


namespace xiaotu {

    
    template <typename Scalar>
    class HermiteDividedDifference
    {
        public:
             /**
             * @brief 构造函数
             * 
             * @param [in] x 采样点 x 列表
             * @param [in] y 对应 x 列表的采样值
             * @param [in] dy 对应 x 列表的一阶导数值
             */
            HermiteDividedDifference(
                    std::vector<Scalar> const & x,
                    std::vector<Scalar> const & y,
                    std::vector<Scalar> const & dy)
            {
                assert(!x.empty());
                assert(x.size() == y.size());
                assert(x.size() == dy.size());
        
                InitTable(x, y, dy);
            }
        
            /**
             * @brief 计算插值
             */
            Scalar Evaluate(Scalar x) const
            {
                int n = mZs.size();
                Scalar result = mTable[0][n - 1]; 
                
                for (int i = n - 2; i >= 0; --i)
                    result = result * (x - mZs[i]) + mTable[0][i];
                return result;
            }

            /**
             * @brief 计算插值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }
            
        private:

            /**
             * @brief 构造广义 Hermite 差商表
             */
            void InitTable(
                std::vector<Scalar> const & x, 
                std::vector<Scalar> const & y,
                std::vector<Scalar> const & dy)
            {
                int num_points = x.size();
                int n = 2 * num_points; // 扩展节点总数
                
                // 1. 构造扩展节点序列 mZs
                mZs.resize(n);
                for (int i = 0; i < num_points; ++i) {
                    mZs[2 * i]     = x[i];
                    mZs[2 * i + 1] = x[i];
                }

                // 2. 初始化差商表矩阵并填充零阶差商（第 0 列为函数值 y）
                mTable.assign(n, std::vector<Scalar>(n, static_cast<Scalar>(0.0)));
                for (int i = 0; i < num_points; ++i) {
                    mTable[2 * i][0]     = y[i];
                    mTable[2 * i + 1][0] = y[i];
                }
        
                // 3. 填充一阶差商（第 1 列）
                // 奇数行(重合点)直接填导数 dy，偶数行用差商公式
                for (int i = 0; i < num_points; ++i) {
                    mTable[2 * i][1] = dy[i];
                    if (i < num_points - 1) {
                        mTable[2 * i + 1][1] = (mTable[2 * i + 2][0] - mTable[2 * i + 1][0])
                                             / (mZs[2 * i + 2] - mZs[2 * i + 1]);
                    }
                }

                // 4. 填充差商表剩余部分
                for (int j = 2; j < n; ++j) {
                    for (int i = 0; i < n - j; ++i) {
                        mTable[i][j] = (mTable[i+1][j-1] - mTable[i][j-1])
                                     / (mZs[i+j] - mZs[i]);
                    }
                }
            }
    
            friend std::ostream & operator << (
                    std::ostream& os,
                    HermiteDividedDifference const & hdd)
            {
                std::ios_base::fmtflags original_flags = os.flags();
                std::streamsize original_precision = os.precision();

                os << std::fixed << std::setprecision(6);
                os << "\n================ Hermite Table ================" << std::endl;

                int n = hdd.mZs.size();
                for (int i = 0; i < n; ++i) {
                    os << "z[" << std::setw(2) << i << "] = " << std::setw(8) << hdd.mZs[i] << " | ";
                    
                    for (int j = 0; j < (n-i); ++j)
                        os << std::setw(12) << hdd.mTable[i][j] << " ";
                    os << std::endl;
                }
                os << "===============================================" << std::endl;

                os.flags(original_flags);
                os.precision(original_precision);
                return os;
            }

        private:
            //! 扩展后的采样点 z (每个原始 x 写两次)
            std::vector<Scalar> mZs;
            //! 广义差商表
            std::vector<std::vector<Scalar>> mTable;
    };

}


#endif
