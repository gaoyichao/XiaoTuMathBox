/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=index
 *
 **************************************************************************** GAO YiChao 2022.0810 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoPoint3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

namespace xiaotu {

    /**
     * @brief 三维射影空间下的点，齐次坐标
     */
    template <typename DataType>
    class HomoPoint3 : public AMatrix<DataType, 4, 1>
    {
        typedef AMatrix<DataType, 4, 1> AVector;
        using AVector::View;

        public:
            HomoPoint3()
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPoint3(DataType _x, DataType _y, DataType _z, DataType _k = 1.0)
            {
                SetValue(_x, _y, _z, _k);
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPoint3(Matrix const & v)
            {
                SetValue(v(0), v(1), v(2), v(3));
            }

            HomoPoint3(HomoPoint3 const & p)
            {
                SetValue(p.x(), p.y(), p.z(), p.k());
            }

            HomoPoint3(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            HomoPoint3 & operator = (HomoPoint3 const & p)
            {
                SetValue(p.x(), p.y(), p.z(), p.k());
                return *this;
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPoint3 & operator = (Matrix const & v)
            {
                SetValue(v(0), v(1), v(2), v(3));
                return *this;
            }

            HomoPoint3 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y, DataType _z, DataType _k)
            {
                x() = _x;
                y() = _y;
                z() = _z;
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
            HomoPoint3 & Normalize()
            {
                *this = Normalization();
                return *this;
            }

            inline bool IsInfinity()
            {
                return 0 == k();
            }

            static HomoPoint3 Infinity()
            {
                return HomoPoint3(1, 0, 0, 0);
            }


        public:
            DataType & x() { return this->At(0); }
            DataType & y() { return this->At(1); }
            DataType & z() { return this->At(2); }
            DataType & k() { return this->At(3); }
            DataType const & x() const { return this->At(0); }
            DataType const & y() const { return this->At(1); }
            DataType const & z() const { return this->At(2); }
            DataType const & k() const { return this->At(3); }
    };
}

