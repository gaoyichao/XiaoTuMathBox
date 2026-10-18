#ifndef XTMB_NUMERICAL_INTEGRAL_HPP
#define XTMB_NUMERICAL_INTEGRAL_HPP


namespace xiaotu {

    /**
     * @brief 梯形法则分段求积
     * 
     * @param [in] y 采样点的被积函数值
     * @param [in] h 采样步长
     */
    template <typename DataType>
    DataType CompositeTrapeRule(std::vector<DataType> const & y, DataType h)
    {
        DataType re = 0;
        size_t n = y.size()-1;

        for (size_t i = 0; i < n; i++) {
            re += y[i] + y[i+1];
        }
        return 0.5 * re * h;
    }

    /**
     * @brief 梯形法则分段求积
     * 
     * @param [in] f 被积函数
     * @param [in] a 积分区间 [a,b]
     * @param [in] b
     * @param [in] n 采样点数量为 n + 1
     */
    template <typename DataType>
    DataType CompositeTrapeRule(std::function<DataType(DataType)> f, DataType a, DataType b, int n)
    {
        std::vector<DataType> ys(n+1);
        DataType step = (b - a) / n;
        for (size_t i = 0; i <= n; i++)
            ys[i] = f(a + i * step);
        return CompositeTrapeRule(ys, step);
    }

    /**
     * @brief Simpson法则分段求积
     * 
     * @param [in] y 采样点的被积函数值, n-1 应当是偶数
     * @param [in] h 采样步长
     */
    template <typename DataType>
    DataType CompositeSimpsonRule(std::vector<DataType> const & y, DataType h)
    {
        DataType re = 0;
        size_t n = y.size()-2;
        assert(1 == (n % 2));

        for (size_t i = 0; i < n; i += 2) {
            re += y[i] + 4 * y[i+1] + y[i+2];
        }
        std::cout << y.back() << std::endl;
        return re * h / 3;
    }
    /**
     * @brief Simpson法则分段求积
     * 
     * @param [in] f 被积函数
     * @param [in] a 积分区间 [a,b]
     * @param [in] b
     * @param [in] n 采样点数量为 n + 1, n 应当是偶数
     */
    template <typename DataType>
    DataType CompositeSimpsonRule(std::function<DataType(DataType)> f, DataType a, DataType b, int n)
    {
        assert(0 == (n % 2));
        std::vector<DataType> ys(n+1);
        DataType step = (b - a) / n;
        for (size_t i = 0; i <= n; i++)
            ys[i] = f(a + i * step);
        return CompositeSimpsonRule(ys, step);
    }


}

#endif