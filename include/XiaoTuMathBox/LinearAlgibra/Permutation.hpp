#ifndef XTMB_LA_PERMUTATION_MATRIX_H
#define XTMB_LA_PERMUTATION_MATRIX_H

#include <cmath>
#include <cassert>
#include <vector>


namespace xiaotu::math {

    class Permutation {
        public:

            Permutation(int n)
            {
                mPerm.resize(n);
                for (int i = 0; i < n; i++)
                    mPerm[i] = i;
            }

            Permutation(std::vector<int> const & perm)
                : mPerm(perm)
            {}


            /**
             * @brief 行交换矩阵, 用于左乘
             */
            template <typename Matrix>
            void LeftMatrix(Matrix & P)
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
            DMatrix<Scalar> LeftMatrix()
            {
                int n = mPerm.size();
                DMatrix<Scalar> re(n, n);
                LeftMatrix(re);
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
             * @brief 列交换矩阵, 用于右乘
             */
            template <typename Matrix>
            void RightMatrix(Matrix & P)
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
             * @brief 行交换矩阵, 用于左乘
             */
            template <typename Scalar>
            DMatrix<Scalar> RightMatrix()
            {
                int n = mPerm.size();
                DMatrix<Scalar> re(n, n);
                RightMatrix(re);
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
