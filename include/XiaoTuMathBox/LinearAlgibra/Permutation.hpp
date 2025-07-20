#ifndef XTMB_LA_PERMUTATION_MATRIX_H
#define XTMB_LA_PERMUTATION_MATRIX_H

#include <cmath>
#include <cassert>
#include <vector>


namespace xiaotu::math {

    /**
     * @brief 置换矩阵
     * 
     * 方阵, 每行/每列只有一个元素为 1，其余全零, 正交
     */
    class Permutation {
        public:

            /**
             * @brief 构造函数
             */
            Permutation()
            {}

            /**
             * @brief 构造函数
             * 
             * @param [in] n 矩阵的维数
             */
            Permutation(int n)
            {
                Resize(n);
            }

            /**
             * @brief 构造函数
             * 
             * @param [in] perm 置换后的行/列映射关系, 由调用者保证元素唯一
             */
            Permutation(std::vector<int> const & perm)
                : mPerm(perm)
            {}


            /**
             * @brief 重置矩阵
             * 
             * @param [in] n 重置后的矩阵维度
             */
            Permutation & Resize(int n)
            {
                mPerm.resize(n);
                for (int i = 0; i < n; i++)
                    mPerm[i] = i;

                return *this;
            }

            /**
             * @brief 置换矩阵的转置
             */
            Permutation Transpose() const
            {
                int n = mPerm.size();

                std::vector<int> tmp;
                tmp.resize(n);

                for (int i = 0; i < n; i++) {
                    tmp[mPerm[i]] = i;
                }
                return Permutation(tmp);
            }

            /**
             * @brief 交换 i,j 行/列
             */
            Permutation & Swap(int i, int j)
            {
                int tmp = mPerm[i];
                mPerm[i] = mPerm[j];
                mPerm[j] = tmp;
                return *this;
            }

            /**
             * @brief 行交换矩阵, 用于左乘
             */
            template <typename Matrix>
            void LeftMatrix(Matrix & P) const
            {
                int n = mPerm.size();
                assert(P.Rows() == P.Cols());
                assert(P.Rows() == n);

                P.Zeroing();
                for (int i = 0; i < n; i++) {
                    P(i, mPerm[i]) = 1;
                }
            }

            /**
             * @brief 行交换矩阵, 用于左乘
             */
            template <typename Scalar>
            DMatrix<Scalar> LeftMatrix() const
            {
                int n = mPerm.size();
                DMatrix<Scalar> re(n, n);
                LeftMatrix(re);
                return re;
            }
 
            /**
             * @brief 列交换矩阵, 用于右乘
             */
            template <typename Matrix>
            void RightMatrix(Matrix & P) const
            {
                int n = mPerm.size();
                assert(P.Rows() == P.Cols());
                assert(P.Rows() == n);

                P.Zeroing();
                for (int i = 0; i < n; i++) {
                    P(mPerm[i], i) = 1;
                }
            }

            /**
             * @brief 列交换矩阵, 用于右乘
             */
            template <typename Scalar>
            DMatrix<Scalar> RightMatrix() const
            {
                int n = mPerm.size();
                DMatrix<Scalar> re(n, n);
                RightMatrix(re);
                return re;
            }

            /**
             * @brief 左乘到矩阵 M 上, M = P * M, 行变换
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void LeftApplyOn(Matrix & M) const
            {
                assert(M.Rows() == mPerm.size());
                DMatrix<typename Matrix::Scalar> tmp = M;
                MatrixRowView rows(tmp, mPerm);
                M = rows;
            }

            /**
             * @brief 左乘到矩阵 A 上
             * 
             * @param [in] A 目标矩阵
             * @return 行变换后的矩阵
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            friend DMatrix<typename Matrix::Scalar>
            operator * (Permutation const & P, Matrix const & A)
            {
                DMatrix<typename Matrix::Scalar> re = A;
                MatrixRowView rows(A, P.mPerm);
                re = rows;
                return re;
            }


            /**
             * @brief 右乘到矩阵 M 上, M = M * P, 列变换
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            void RightApplyOn(Matrix & M) const
            {
                assert(M.Cols() == mPerm.size());
                DMatrix<typename Matrix::Scalar> tmp = M;
                MatrixColView cols(tmp, mPerm);
                M = cols;
            }

            /**
             * @brief 右乘到矩阵 A 上
             * 
             * @param [in] A 目标矩阵
             * @return 列变换后的矩阵
             */
            template <typename Matrix, bool MIsMatrix = Matrix::IsMatrix>
            friend DMatrix<typename Matrix::Scalar>
            operator * (Matrix const & A, Permutation const & P)
            {
                DMatrix<typename Matrix::Scalar> re = A;
                MatrixColView cols(A, P.mPerm);
                re = cols;
                return re;
            }

        private:
            std::vector<int> mPerm;
    };

}


#endif
