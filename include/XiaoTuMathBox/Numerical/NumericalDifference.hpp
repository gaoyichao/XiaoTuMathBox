#ifndef XTMB_NUMERICAL_DIFFERENCE_HPP
#define XTMB_NUMERICAL_DIFFERENCE_HPP


namespace xiaotu {

    /**
     * @brief 前向差分
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType ForwardDifference(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (f(x0 + h) - f(x0)) / h;
    }

    /**
     * @brief 后向差分
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType BackwardDifference(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (f(x0 - h) - f(x0)) / (-h);
    }

    /**
     * @brief 三点中心公式
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType MidThreePoints(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (f(x0 + h) - f(x0 - h)) / h * 0.5;
    }


    /**
     * @brief 三点端点公式
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType EndThreePoints(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (-3 * f(x0) + 4 * f(x0 + h) - f(x0 + 2 * h)) / h * 0.5;
    }


    /**
     * @brief 五点中心公式
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType MidFivePoints(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (f(x0 - 2 * h) - 8 *f(x0 - h) + 8 * f(x0 + h) - f(x0 + 2*h)) / h / 12;
    }

    /**
     * @brief 五点端点公式
     * 
     * @param [in] f 目标函数
     * @param [in] x0 考察点
     * @param [in] h 差分步长
     */
    template <typename DataType>
    DataType EndFivePoints(std::function<DataType(DataType)> f, DataType x0, DataType h)
    {
        return (-25 * f(x0) + 48 * f(x0+h) - 36*f(x0+2*h) + 16*f(x0+3*h) - 3*f(x0+4*h)) / h / 12;
    }


}

#endif

