/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=index
 *
 **************************************************************************** GAO YiChao 2022.0810 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoPlane3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

namespace xiaotu {

    /**
     * @brief 三维射影空间下的平面，齐次坐标
     */
    template <typename DataType>
    class HomoPlane3 : public AMatrix<DataType, 4, 1>
    {
        typedef AMatrix<DataType, 4, 1> AVector;
        using AVector::View;

        public:
            HomoPlane3()
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPlane3(DataType _a, DataType _b, DataType _c, DataType _d = 1.0)
            {
                SetValue(_a, _b, _c, _d);
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPlane3(Matrix const & v)
            {
                SetValue(v(0), v(1), v(2), v(3));
            }

            HomoPlane3(HomoPlane3 const & p)
            {
                SetValue(p.a(), p.b(), p.c(), p.d());
            }

            HomoPlane3(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            HomoPlane3 & operator = (HomoPlane3 const & p)
            {
                SetValue(p.a(), p.b(), p.c(), p.d());
                return *this;
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoPlane3 & operator = (Matrix const & v)
            {
                SetValue(v(0), v(1), v(2), v(3));
                return *this;
            }

            HomoPlane3 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c, DataType _d)
            {
                a() = _a;
                b() = _b;
                c() = _c;
                d() = _d;
            }

            //! @brief 归一化，不修改对象本身
            //!
            //! @param [in] tolerance 分母为0的判定容忍度
            AVector Normalization(DataType tolerance = SMALL_VALUE) const
            {
                DataType norm = this->Norm();
                if (norm < tolerance) {
                    return {0, 0, 0, 1};
                } else {
                    DataType norm_inv = 1.0 / norm;
                    return (*this) * norm_inv;
                }
            }
            //! @brief 归一化，修改对象本身
            //!
            //! @param [in] tolerance 分母为0的判定容忍度
            HomoPlane3 & Normalize(DataType tolerance = SMALL_VALUE)
            {
                *this = Normalization(tolerance);
                return *this;
            }

            inline bool IsInfinity()
            {
                return d() != 0 && a() == 0 && b() == 0 && c() == 0;
            }

            static HomoPlane3 Infinity()
            {
                return HomoPlane3(0, 0, 0, 1);
            }

        public:
            DataType & a() { return this->At(0); }
            DataType & b() { return this->At(1); }
            DataType & c() { return this->At(2); }
            DataType & d() { return this->At(3); }
            DataType const & a() const { return this->At(0); }
            DataType const & b() const { return this->At(1); }
            DataType const & c() const { return this->At(2); }
            DataType const & d() const { return this->At(3); }
    };

}

