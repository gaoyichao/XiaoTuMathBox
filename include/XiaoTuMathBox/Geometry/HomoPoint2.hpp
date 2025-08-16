/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoPoint2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

namespace xiaotu::math {

    /**
     * @brief 二维摄影空间下的向量
     */
    template <typename DataType>
    class HomoPoint2 : public AMatrix<DataType, 3, 1>
    {
        typedef AMatrix<DataType, 3, 1> AVector;
        using AVector::View;
        
        public:
            HomoPoint2()
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoPoint2(DataType _x, DataType _y, DataType _k = 1.0)
            {
                SetValue(_x, _y, _k);
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPoint2(Matrix const & v)
            {
                SetValue(v(0), v(1), v(2));
            }

            HomoPoint2(HomoPoint2 const & p)
            {
                SetValue(p.x(), p.y(), p.k());
            }

            HomoPoint2(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            HomoPoint2 & operator = (HomoPoint2 const & p)
            {
                SetValue(p.x(), p.y(), p.k());
                return *this;
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPoint2 & operator = (Matrix const & v)
            {
                SetValue(v(0), v(1), v(2));
                return *this;
            }

            HomoPoint2 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y, DataType _k)
            {
                x() = _x;
                y() = _y;
                k() = _k;
            }

            /**
             * @brief 归一化, 构造新对象，不修改对象本身
             */
            AVector Normalization() const
            {
                if (0 == k())
                    return *this;

                DataType k_inv = 1.0 / k();
                return *this * k_inv;
            }

            /**
             * @brief 归一化，修改对象本身
             */
            HomoPoint2 & Normlize()
            {
                *this = Normalization();
                return *this;
            }

            inline bool IsInfinity()
            {
                return 0 == k();
            }

            static HomoPoint2 Infinity()
            {
                return HomoPoint2(1, 0, 0);
            }

        public:
            inline DataType & x() { return this->At(0); }
            inline DataType & y() { return this->At(1); }
            inline DataType & k() { return this->At(2); }
            inline DataType const & x() const { return this->At(0); }
            inline DataType const & y() const { return this->At(1); }
            inline DataType const & k() const { return this->At(2); }
    };

}

