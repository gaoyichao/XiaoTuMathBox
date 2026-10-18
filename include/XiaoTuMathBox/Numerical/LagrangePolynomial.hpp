/************************************************************************
 * 
 * 拉格朗日多项式(Lagrange Polynomial)
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=拉格朗日多项式
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_H
#define XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_H

namespace xiaotu {

    /**
     * @brief 拉格朗日基函数
     * 
     * L_{n,i}(x) = \frac{(  x - x_0)\cdots(  x - x_{i-1})(  x - x_{i+1}) \cdots (  x - x_n)}
     *                   {(x_i - x_0)\cdots(x_i - x_{i-1})(x_i - x_{i+1}) \cdots (x_i - x_n)}
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] i 参考采样点索引
     * @param [in] x 目标插值点
     * @return 拉格朗日基函数值
     */
    template <typename DataType>
    DataType LagrangeBasis(std::vector<DataType> const & x_nodes,
                           size_t i, DataType x)
    {
        size_t n = x_nodes.size();
        DataType Li = 1;
        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;
            Li *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
        }
        return Li;
    }

    /**
     * @brief 拉格朗日基函数
     * 
     * L_{n,i}(x) = \frac{(  x - x_0)\cdots(  x - x_{i-1})(  x - x_{i+1}) \cdots (  x - x_n)}
     *                   {(x_i - x_0)\cdots(x_i - x_{i-1})(x_i - x_{i+1}) \cdots (x_i - x_n)}
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] i 参考采样点索引
     * @return 拉格朗日基多项式
     */
    template <typename DataType>
    Polynomial<DataType> LagrangeBasis(std::vector<DataType> const & x_nodes, size_t i)
    {
        size_t n = x_nodes.size();
        std::vector<DataType> L_i = {1.0}; 
        // L_i(x) = \prod_{j \neq i} (x - x_j) / (x_i - x_j)
        DataType factor = 1.0;
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
            
            factor *= (x_nodes[i] - x_nodes[j]);
        }
        
        Polynomial<DataType> re(L_i);
        factor = Traits<Polynomial<DataType>>::one / factor;
        for (int j = 0; j < n; ++j)
            re[j] *= factor;
        return re;
    } 

    /**
     * @brief 拉格朗日基函数在 x_i 处的导数
     * 
     * L_{n,i}/(x) = \sum_{i \neq j} frac{1}{x_i - x_j}
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] i 参考采样点索引
     * @return 拉格朗日基函数在 x_i 处的一阶导数值
     */
    template <typename DataType>
    DataType LagrangeBasisDerivative(std::vector<DataType> const & x_nodes,
                           size_t i)
    {
        size_t n = x_nodes.size();
        DataType dLi = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;
            dLi += 1 / (x_nodes[i] - x_nodes[j]);
        }
        return dLi;
    }

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
            result += y_nodes[i] * LagrangeBasis(x_nodes, i, x);
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
    template <typename DataType>
    Polynomial<DataType> LagrangePolynomial(std::vector<DataType> const & x_nodes, 
                                            std::vector<DataType> const & y_nodes)
    {
        size_t n = x_nodes.size();
        assert(y_nodes.size() == n);

        Polynomial<DataType> re;
        re.ReAlloc(n-1);

        for (int i = 0; i < n; ++i) {
            Polynomial<DataType> L_i = LagrangeBasis(x_nodes, i);
            for (int k = 0; k < n; ++k)
                re[k] += y_nodes[i] * L_i[k];
        }

        return re;
    }


    /**
     * @brief Hermite 插值多项式的朴素实现
     * 
     * @param [in] x_nodes 采样点 x 坐标
     * @param [in] y_nodes 采样点 y 坐标
     * @param [in] d_nodes 采样点的一阶导数
     * @param [in] x 插值点
     */
    template <typename DataType>
    DataType HermiteInterpolation(std::vector<DataType> const & x_nodes, 
                                  std::vector<DataType> const & y_nodes, 
                                  std::vector<DataType> const & d_nodes, 
                                  DataType x)
    {
        size_t n = x_nodes.size();
        assert(y_nodes.size() == n);
        assert(d_nodes.size() == n);

        DataType re = 0;
        for (size_t i = 0; i <= n; ++i) {
            DataType L_i = LagrangeBasis(x_nodes, i, x);
            DataType dL_i = LagrangeBasisDerivative(x_nodes, i);

            // H_{n,i}(x) = [1 - 2 * L'_{n,i}(x_i) * (x - x_i)] * L^2_{n,i}(x)
            DataType H = (1 - 2 * dL_i * (x - x_nodes[i])) * (L_i * L_i);

            // \hat{H}_{n,i}(x) = (x - x_i) * L^2_{n,i}(x)
            DataType hat_H = (x - x_nodes[i]) * (L_i * L_i);

            re += y_nodes[i] * H + d_nodes[i] * hat_H;
        }

        return re;
    }



}


#endif