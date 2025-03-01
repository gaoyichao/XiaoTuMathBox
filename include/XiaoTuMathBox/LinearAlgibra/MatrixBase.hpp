#ifndef XTMB_LA_MATRIX_BASE_H
#define XTMB_LA_MATRIX_BASE_H

#include <cmath>
#include <cassert>
#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename Derived> 
    class MatrixBase
    {
        public:
            typedef typename Traits<Derived>::Scalar Scalar;
            constexpr static bool IsMatrix = true;
            constexpr static EAlignType Align = Traits<Derived>::Align;
            constexpr static EStoreType Store = Traits<Derived>::Store;
        public:

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            MatrixBase & operator = (std::initializer_list<Scalar> && li)
            {
                assert(li.size() == NumDatas());
                Assign(0, 0, Rows(), Cols(), li.begin());
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            //!
            //! 调用该函数之前，需要保证行数和列数能够容纳下 li 的所有数据。
            //! 不处理 li 没有覆盖到的部分。
            //!
            //! @param [in] li 用于拷贝的初始化列表
            MatrixBase & operator = (std::initializer_list<std::initializer_list<Scalar>> && li)
            {
                Assign(std::move(li));
                return *this;
            }

            //! @brief 拷贝赋值
            template <typename Mat>
            MatrixBase & operator = (Mat const & mv)
            {
                Assign(mv);
                return *this;
            }

            //! @brief 深度拷贝 m
            //!
            //! @param [in] m 将矩阵 m 中的值拷贝过来
            template <typename Mat>
            void Assign(Mat const & m)
            {
                assert(Rows() == m.Rows() && Cols() == m.Cols());

                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = m(ridx, cidx);
            }

            //! @brief 拷贝矩阵
            //!
            //! @param [in] r 子阵起始行索引
            //! @param [in] c 子阵起始列索引
            template <typename Mat>
            void Assign(int r, int c, Mat const & m)
            {
                assert(r >= 0 && r < Rows());
                assert(c >= 0 && c < Cols());

                for (int ridx = 0; ridx < m.Rows(); ridx++) {
                    int i = r + ridx;
                    if (i > Rows())
                        break;
                    for (int cidx = c; cidx < m.Cols(); cidx++) {
                        int j = c + cidx;
                        if (j > Cols())
                            break;
                        At(i, j) = m(ridx, cidx);
                    }
                }
            }

            template <typename Vec>
            void AssignCol(int c, Vec const & m)
            {
                assert(c >= 0 && c < Cols());
                int rows = Rows() < m.NumDatas() ? Rows() : m.NumDatas();

                for (int ridx = 0; ridx < rows; ridx++) {
                    At(ridx, c) = m(ridx);
                }
            }

            template <typename Vec>
            void AssignRow(int r, Vec const & m)
            {
                assert(r >= 0 && r < Rows());
                int cols = Cols() < m.NumDatas() ? Cols() : m.NumDatas();

                for (int cidx = 0; cidx < cols; cidx++) {
                    At(cidx, r) = m(cidx);
                }
            }

            //! @brief 深度拷贝 vlist 到 s 所指的连续内存中
            //!
            //! @param [in] s 目标起始索引
            //! @param [in] n 拷贝数据长度
            //! @param [in] vlist 用户维护的内存
            void Assign(int s, int n, Scalar const * vlist)
            {
                assert((s + n) <= NumDatas());
                Scalar * start = Ptr(s);
                for (int i = 0; i < n; i++)
                    start[i] = vlist[i];
            }

            //! @brief 深度拷贝 vlist 到 ridx, cidx 所指的子阵中
            //!
            //! @param [in] r 子阵起始行索引
            //! @param [in] c 子阵起始列索引
            //! @param [in] m 子阵行数
            //! @param [in] n 子阵列数
            //! @param [in] vlist 用户维护的内存
            void Assign(int r, int c, int m, int n, Scalar const * vlist)
            {
                int rend = r + m;
                int cend = c + n;

                assert(r >= 0 && rend <= Rows());
                assert(c >= 0 && cend <= Cols());

                int vidx = 0;
                for (int ridx = r; ridx < rend; ridx++)
                    for (int cidx = c; cidx < cend; cidx++)
                        At(ridx, cidx) = vlist[vidx++];
            }

            //! @brief 深度拷贝 li
            //!
            //! 调用该函数之前，需要保证行数和列数能够容纳下 li 的所有数据。
            //! 不处理 li 没有覆盖到的部分。
            //!
            //! @param [in] li 用于拷贝的初始化列表
            void Assign(std::initializer_list<std::initializer_list<Scalar>> && li)
            {
                assert(Rows() >= li.size());

                int ridx = 0;
                for (auto rit = li.begin(); rit != li.end(); rit++, ridx++) {
                    assert(Cols() >= rit->size());
                    int cidx = 0;
                    for (auto cit = rit->begin(); cit != rit->end(); cit++, cidx++)
                        At(ridx, cidx) = *cit;
                }
            }

            //! @brief 填充整个矩阵
            //!
            //! @param [in] 填充数值
            void Full(Scalar const & v)
            {
                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = v;
            }

            //! @brief 将矩阵中所有元素都置为零
            void Zeroing()
            {
                for (int i = 0; i < Rows(); i++)
                    for (int j = 0; j < Cols(); j++)
                        At(i, j) = 0;
            }

            //! @brief 单位化, 将矩阵改写成单位矩阵
            void Identity()
            {
                Zeroing();
                int N = Rows() < Cols() ? Rows() : Cols();
                for (int i = 0; i < N; ++i) {
                    At(i, i) = 1;
                }
            }

            //! @brief 交换 i, j 两行
            void RowSwap(int i, int j)
            {
                assert(i != j);
                assert(i < Rows() && j < Rows());
                for (int k = 0; k < Cols(); k++)
                    std::swap(At(i, k), At(j, k));
            }

            //! @brief 交换 i, j 两列
            void ColSwap(int i, int j)
            {
                assert(i != j);
                assert(i < Cols() && j < Cols());
                for (int k = 0; k < Rows(); k++)
                    std::swap(At(k, i), At(k, j));
            }

            MatrixBase & Normalize()
            {
                Scalar norm = this->Norm();
                if (std::abs(norm) < SMALL_VALUE) {
                    for (int i = 0; i < Rows(); i++)
                        for (int j = 0; j < Cols(); j++)
                            At(i, j) = 0;
                } else {
                    Scalar norm_inv = 1.0 / norm;
                    for (int i = 0; i < Rows(); i++)
                        for (int j = 0; j < Cols(); j++)
                            At(i, j) *= norm_inv;
                }

                return *this;
            }


        public:
            ////////////////////////////////////////////////////////
            //
            //  一些聚合运算
            //
            ////////////////////////////////////////////////////////

            //! @brief 所有元素的平方和
            Scalar SquaredNorm() const
            {
                return Dot(*this);
            }

            //! @brief 二范数, 模
            Scalar Norm() const
            {
                return std::sqrt(SquaredNorm());
            }

        public:
            ////////////////////////////////////////////////////////
            //
            //  与其他矩阵的运算
            //
            ////////////////////////////////////////////////////////

            //! @brief 点乘
            template <typename MType>
            Scalar Dot(MType const & m) const
            {
                assert(NumDatas() == m.NumDatas());
                Scalar re = 0;
                for (int i = 0; i < NumDatas(); i++)
                    re += At(i) * m(i);
                return re;
            }

        public:
            ////////////////////////////////////////////////////////
            //
            //  访问矩阵元素
            //
            ////////////////////////////////////////////////////////

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                return (EAlignType::eRowMajor == Align)
                      ? Cols() * row + col
                      : Rows() * col + row;
            }

            //! @brief 获取指定位置的元素指针
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素指针
            inline Scalar * Ptr(int idx = 0)
            {
                return derived().StorBegin() + idx;
            }

            //! @brief 获取指定位置的元素指针
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素指针
            inline Scalar const * Ptr(int idx = 0) const
            {
                return derived().StorBegin() + idx;
            }

            //! @brief 获取指定位置的元素指针
            //!
            //! @param [in] row 元素的行索引
            //! @param [in] col 元素的列索引
            //! @return 元素指针
            inline Scalar * Ptr(int row, int col)
            {
                return derived().StorBegin() + Idx(row, col);
            }

            //! @brief 获取指定位置的元素指针
            //!
            //! @param [in] row 元素的行索引
            //! @param [in] col 元素的列索引
            //! @return 元素指针
            inline Scalar const * Ptr(int row, int col) const
            {
                return derived().StorBegin() + Idx(row, col);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar & At(int idx)
            {
                return *(Ptr(idx));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar const & At(int idx) const
            {
                return *(Ptr(idx));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar & At(int row, int col)
            {
                return *(Ptr(row, col));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar const & At(int row, int col) const
            {
                return *(Ptr(row, col));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar & operator() (int idx)
            {
                return At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar const & operator() (int idx) const
            {
                return At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar & operator() (int row, int col)
            {
                return At(row, col);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar const & operator() (int row, int col) const
            {
                return At(row, col);
            }

        public:
            //! @brief 获取矩阵行数
            inline int Rows() const { return derived().Rows(); }
            //! @brief 获取矩阵列数
            inline int Cols() const { return derived().Cols(); }
            //! @brief 获取矩阵元素数量
            inline int NumDatas() const { return Rows() * Cols(); }
            //! @brief 获取矩阵总字节数
            inline int NumBytes() const { return NumDatas() * sizeof(Scalar); }
            //! @brief 获取矩阵元素字节数
            inline int DataBytes() const { return sizeof(Scalar); }
            //! @brief 获取矩阵的存储方式
            inline EAlignType AlignType() const { return Align; }
 
            friend std::ostream & operator << (std::ostream & s, MatrixBase const & m)
            {
                s << std::endl;
                for (int ridx = 0; ridx < m.Rows(); ridx++) {
                    s << m(ridx, 0);
                    for (int cidx = 1; cidx < m.Cols(); cidx++)
                        s << "\t" << m(ridx, cidx);
                    s << std::endl;
                }
                return s;
            }
       
        private:
            Derived & derived() { return *static_cast<Derived*>(this); }
            Derived const & derived() const { return *static_cast<const Derived*>(this); }
    };

}

#endif
