#ifndef XTMB_LA_GIVENS_MATRIX_H
#define XTMB_LA_GIVENS_MATRIX_H

#include <cmath>
#include <cassert>
#include <vector>

namespace xiaotu::math {

    template <typename Scalar>
    class Givens {
        public:

            /**
             * @brief 拷贝构造
             */
            Givens(Givens const & other)
                : i_(other.i_), j_(other.j_), c_(other.c_), s_(other.s_)
            {}

            Givens & operator = (Givens const & other)
            {
                i_ = other.i_;
                j_ = other.j_;
                c_ = other.c_;
                s_ = other.s_;
                return *this;
            }

            /**
             * @brief 构造函数
             *
             * @param [in] i 旋转 i_-j_ 平面
             * @param [in] j 旋转 i_-j_ 平面
             * @param [in] theta 旋转弧度
             */
            Givens(int i, int j, Scalar theta)
                : i_(i), j_(j)
            {
                c_ = std::cos(theta);
                s_ = std::sin(theta);
            }

            /**
             * @brief 构造函数
             * 
             * @param [in] i 旋转 i_-j_ 平面
             * @param [in] j 旋转 i_-j_ 平面
             * @param [in] vi 向量 v 中第 i_ 行元素值
             * @param [in] vj 向量 v 中第 j_ 行元素值
             */
            Givens(int i, int j, Scalar vi, Scalar vj)
            {
                Set(i, j, vi, vj);
            }


            /**
             * @brief 构建 Givens 矩阵 G(i_,j_, \theta)
             *
             * @param [in] i_ 旋转 i_-j_ 平面
             * @param [in] j_ 旋转 i_-j_ 平面
             * @param [in] v 向量 v
             */
            template <typename Vector>
            Givens(int i, int j, Vector const & v)
                : i_(i), j_(j)
            {
                assert(v.NumDatas() >= 2);

                Scalar vi = v(i_);
                Scalar vj = v(j_);

                Set(i, j, vi, vj);
            }

            /**
             * @brief 重置 Givens(i_, j_, \theta) 矩阵参数 
             * 
             * @param [in] i 旋转 i_-j_ 平面
             * @param [in] j 旋转 i_-j_ 平面
             * @param [in] vi 向量 v 中第 i_ 行元素值
             * @param [in] vj 向量 v 中第 j_ 行元素值
             */
            void Set(int i, int j, Scalar vi, Scalar vj)
            {
                i_ = i;
                j_ = j;

                if (std::abs(vj) < SMALL_VALUE) {
                    c_ = std::copysign(1, vi);
                    s_ = 0;
                } else {
                    if (std::abs(vj) > std::abs(vi)) {
                        Scalar tau = vi / vj;
                        s_ = std::copysign(1 / std::sqrt(1 + tau * tau), vj);
                        c_ = s_ * tau;
                    } else {
                        Scalar tau = vj / vi;
                        c_ = std::copysign(1 / std::sqrt(1 + tau * tau), vi);
                        s_ = c_ * tau;
                    }
                }
            }

            /**
             * @brief 转换成矩阵的形式
             */
            template <typename Matrix>
            void ToMatrix(Matrix & G)
            {
                assert(G.Rows() == G.Cols());
                assert(i_ < G.Rows() && j_ < G.Rows() && G.Rows() >= 2);

                G.Identity();
                G(i_,i_) = c_; G(i_,j_) = s_;
                G(j_,i_) = -s_; G(j_,j_) = c_;
            }
    
            /**
             * @brief 转换成 n x n 的矩阵
             * 
             * @param [in] n 输出矩阵大小
             */
            DMatrix<Scalar> ToMatrix(int n)
            {
                assert(i_ < n && j_ < n && n >= 2);
                DMatrix<Scalar> re(n, n);
                ToMatrix(re);
                return re;
            }

            /**
             * @brief 作用到矩阵 M 上, M = G * M
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void LeftApplyOn(Matrix & M) const
            {
                int rows = M.Rows();
                int cols = M.Cols();
                assert(i_ < rows && j_ < rows && rows >= 2);

                for (int cidx = 0; cidx < cols; cidx++) {
                    Scalar vi =  M(i_, cidx) * c_ + M(j_, cidx) * s_;
                    Scalar vj = -M(i_, cidx) * s_ + M(j_, cidx) * c_;

                    M(i_, cidx) = vi;
                    M(j_, cidx) = vj;
                }
            }

            /**
             * @brief 作用到矩阵 M 上, M = G^T * M
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void TLeftApplyOn(Matrix & M) const
            {
                int rows = M.Rows();
                int cols = M.Cols();
                assert(i_ < rows && j_ < rows && rows >= 2);

                for (int cidx = 0; cidx < cols; cidx++) {
                    Scalar vi = M(i_, cidx) * c_ - M(j_, cidx) * s_;
                    Scalar vj = M(i_, cidx) * s_ + M(j_, cidx) * c_;

                    M(i_, cidx) = vi;
                    M(j_, cidx) = vj;
                }
            }

            /**
             * @brief 作用到矩阵 M 上, M = M * G
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void RightApplyOn(Matrix & M) const
            {
                int rows = M.Rows();
                int cols = M.Cols();
                assert(i_ < cols && j_ < cols && cols >= 2);

                for (int ridx = 0; ridx < rows; ridx++) {
                    Scalar vi = M(ridx, i_) * c_ - M(ridx, j_) * s_;
                    Scalar vj = M(ridx, i_) * s_ + M(ridx, j_) * c_;

                    M(ridx, i_) = vi;
                    M(ridx, j_) = vj;
                }
            }
            /**
             * @brief 作用到矩阵 M 上, M = M * G^T
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void TRightApplyOn(Matrix & M) const
            {
                int rows = M.Rows();
                int cols = M.Cols();
                assert(i_ < cols && j_ < cols && cols >= 2);

                for (int ridx = 0; ridx < rows; ridx++) {
                    Scalar vi =  M(ridx, i_) * c_ + M(ridx, j_) * s_;
                    Scalar vj = -M(ridx, i_) * s_ + M(ridx, j_) * c_;

                    M(ridx, i_) = vi;
                    M(ridx, j_) = vj;
                }
            }


            /**
             * @brief 左乘到矩阵 A 上
             * 
             * @param [in] A 目标矩阵
             * @return Givens 变换后的矩阵
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            friend DMatrix<Scalar> operator * (Givens const & G, Matrix const & A)
            {
                DMatrix<Scalar> re = A;
                G.LeftApplyOn(re);
                return re;
            }

            /**
             * @brief 右乘到矩阵 A 上
             * 
             * @param [in] A 目标矩阵
             * @return Givens 变换后的矩阵
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            friend DMatrix<Scalar> operator * (Matrix const & A, Givens const & G)
            {
                DMatrix<Scalar> re = A;
                G.RightApplyOn(re);
                return re;
            }

            Scalar c() const { return c_; }
            Scalar s() const { return s_; }
        private:
            //! @brief 旋转 i_-j_ 平面
            int i_{0};
            //! @brief 旋转 i_-j_ 平面
            int j_{1};
            //! @brief cos(\theta)
            Scalar c_{1};
            //! @brief sin(\theta)
            Scalar s_{0};
    };


    //! @brief 构建 Givens 矩阵 G(i,j, \theta)
    //!
    //! @param [in] i 旋转 i-j 平面
    //! @param [in] j 旋转 i-j 平面
    //! @param [in] v 向量 v
    //! @return G 输出的 Givens 矩阵, 使得 Gv(j) = 0
    template <typename Vector>
    DMatrix<typename Vector::Scalar>
    GivensMatrix(int i, int j, Vector const & v)
    {
        using Scalar = typename Vector::Scalar;
        assert(v.Rows() >= 2 && v.Cols() == 1);

        Givens<Scalar> g(i, j, v);
        return g.ToMatrix(v.Rows());
    }

}

#endif
