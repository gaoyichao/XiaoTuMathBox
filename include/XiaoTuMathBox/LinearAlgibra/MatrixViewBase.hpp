#ifndef XTMB_LA_MATRIX_VIEW_INTERFACE_H
#define XTMB_LA_MATRIX_VIEW_INTERFACE_H

#include <stdint.h>
#include <sys/types.h>
#include <cassert>
#include <algorithm>

namespace xiaotu {
namespace math {

    //! @brief 矩阵视图的基类
    template <typename Derived>
    class MatrixViewBase {
        public:
            typedef typename traits<Derived>::Scalar Scalar;

        public:
            virtual ~MatrixViewBase() {}
        public:
            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            virtual int Idx(int row, int col) const = 0;
            virtual Scalar * StorBegin() = 0;
            virtual Scalar const * StorBegin() const = 0;

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar * Ptr(int idx = 0)
            {
                return StorBegin() + idx;
            }

            inline Scalar const * Ptr(int idx = 0) const
            {
                return StorBegin() + idx;
            }

            inline Scalar * Ptr(int row, int col)
            {
                return StorBegin() + Idx(row, col);
            }

            inline Scalar const * Ptr(int row, int col) const
            {
                return StorBegin() + Idx(row, col);
            }

        public:

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
            //! @brief 填充整个矩阵
            //!
            //! @param [in] 填充数值
            void Full(Scalar const & v)
            {
                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = v;
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
            void Assign(std::initializer_list<std::initializer_list<Scalar>> li)
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

        public:
            //! @brief 交换 i, j 两行
            inline void RowSwap(int i, int j)
            {
                assert(i != j);
                assert(i < Rows() && j < Rows());
                for (int k = 0; k < Cols(); k++)
                    std::swap(At(i, k), At(j, k));
            }

            //! @brief 交换 i, j 两列
            inline void ColSwap(int i, int j)
            {
                assert(i != j);
                assert(i < Cols() && j < Cols());
                for (int k = 0; k < Rows(); k++)
                    std::swap(At(k, i), At(k, j));
            }

            friend std::ostream & operator << (std::ostream & s, Derived const & m)
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

        public:
            //! @brief 获取矩阵行数
            virtual int Rows() const = 0;
            //! @brief 获取矩阵列数
            virtual int Cols() const = 0;
            //! @brief 获取矩阵元素数量
            inline int NumDatas() const { return Rows() * Cols(); }
            //! @brief 获取矩阵总字节数
            inline int NumBytes() const { return NumDatas() * sizeof(Scalar); }
            //! @brief 获取矩阵元素字节数
            inline int DataBytes() const { return sizeof(Scalar); }
    };
}
}

#endif
