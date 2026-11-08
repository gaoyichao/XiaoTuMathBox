#ifndef XTMB_NUMERICAL_LINEAR_SYSTEM
#define XTMB_NUMERICAL_LINEAR_SYSTEM

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu {


    template <typename Scalar>
    class LinearSystem
    {
        public:
            template <typename XMat, typename AMat, typename BMat,
                      bool XIsMatrix = XMat::IsMatrix,
                      bool AIsMatrix = AMat::IsMatrix,
                      bool BIsMatrix = BMat::IsMatrix>
            LinearSystem(XMat const & x0, AMat const & A, BMat const & B)
                : _x(x0), _A(A), _B(B)
            {
                assert(_x.Cols() == 1);
                assert(_x.Rows() == _A.Cols() && _A.Cols() == _A.Rows());
                assert(_B.Rows() == _x.Rows());
            }


            /**
             * @brief 微分方程 dx = Ax + Bu
             * 
             * @param [in] pu 控制量指针, nullptr 表示没有控制量
             */
            template <typename UMat, bool UIsMatrix = UMat::IsMatrix>
            DMatrix<Scalar> DotX(UMat const * pu = nullptr)
            {
                DMatrix<Scalar> re = _A * _x;
                if (nullptr != pu) {
                    assert(pu->Rows() == _B.Cols());
                    re += _B * (*pu);
                }
                return re;
            }

            /**
             * @brief 按照欧拉法推进一个步长 dt
             * 
             * @param [in] dt 仿真步长
             * @param [in] pu 控制量指针, nullptr 表示没有控制量
             */
            template <typename UMat, bool UIsMatrix = UMat::IsMatrix>
            DMatrix<Scalar> EulerStep(double dt, UMat const * pu = nullptr)
            {
                auto dx = DotX(pu);
                return _x + dx * dt;
            }

            /**
             * @brief 更新系统状态
             */
            template <typename XMat, bool XIsMatrix = XMat::IsMatrix>
            LinearSystem & UpdateState(XMat const & x)
            {
                assert(x.Cols() == _x.Cols() && x.Rows() == _x.Rows());
                _x = x;
                return *this;
            }

            /**
             * @brief 当前系统状态
             */
            DMatrix<Scalar> const & State() const
            {
                return _x;
            }

        private:
            DMatrix<Scalar> _x;
            DMatrix<Scalar> _A;
            DMatrix<Scalar> _B;

    };

}

#endif

