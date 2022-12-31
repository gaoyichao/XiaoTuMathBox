#ifndef XTMB_LA_MATRIX_VIEW_H
#define XTMB_LA_MATRIX_VIEW_H

#include <iostream>
#include <XiaoTuMathBox/LinearAlgibra/MatrixViewBase.hpp>

namespace xiaotu {
namespace math {

    template <typename T>
    class MatrixView : public MatrixViewBase
    {
        public:
            //! @brief 构造函数
            //!
            //! @param [in] buffer 数据缓存
            //! @param [in] num 数据长度
            MatrixView(T * buffer, int num)
                : MatrixViewBase((uint8_t*)buffer, num * sizeof(T))
            {
                mDataBytes = sizeof(T);
            }

            //! @brief 浅拷贝构造
            MatrixView(MatrixView const & mv)
                : MatrixViewBase(mv)
            {
                mDataBytes = sizeof(T);
            }

            //! @brief 浅拷贝赋值
            MatrixView & operator = (MatrixView const & mv)
            {
                mStorBegin = mv.mStorBegin;
                mBytes = mv.mBytes;
                Reshape(mv.mNumRows, mv.mNumCols, mv.mDataBytes, mv.mOption);
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            MatrixView & operator = (std::initializer_list<T> li)
            {
                assert(li.size() == NumDatas());
                Assign(0, 0, mNumRows, mNumCols, li.begin());
                return *this;
            }

            //! @brief 重置矩阵形状，不重新分配内存
            //!
            //! @param [in] rows 新形状的行数
            //! @param [in] cols 新形状的列数
            //! @param [in] option 缓存的解析方式，默认列优先解析
            inline void Reshape(int rows, int cols, EStorageOptions option = eColMajor)
            {
                MatrixViewBase::Reshape(rows, cols, sizeof(T), option);
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素地址
            inline T * Ptr(int idx)
            {
                return MatrixViewBase::Ptr<T>(idx);
            }

            inline T const * Ptr(int idx) const
            {
                return MatrixViewBase::Ptr<T>(idx);
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素地址
            inline T * Ptr(int row, int col)
            {
                return MatrixViewBase::Ptr<T>(row, col);
            }

            inline T const * Ptr(int row, int col) const
            {
                return MatrixViewBase::Ptr<T>(row, col);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline T & At(int idx)
            {
                return *(Ptr(idx));
            }

            inline T const & At(int idx) const
            {
                return *(Ptr(idx));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline T & At(int row, int col)
            {
                return *(Ptr(row, col));
            }

            inline T const & At(int row, int col) const
            {
                return *(Ptr(row, col));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline T & operator() (int idx)
            {
                return At(idx);
            }

            inline T const & operator() (int idx) const
            {
                return At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline T & operator() (int row, int col)
            {
                return At(row, col);
            }

            inline T const & operator() (int row, int col) const
            {
                return At(row, col);
            }

        public:
            //! @brief 填充整个矩阵
            //!
            //! @param [in] 填充数值
            inline void Full(T const & v)
            {
                MatrixViewBase::Full<T>(v);
            }

            //! @brief 深度拷贝 vlist 到 s 所指的连续内存中
            //!
            //! @param [in] s 目标起始索引
            //! @param [in] n 拷贝数据长度
            //! @param [in] vlist 用户维护的内存
            void Assign(int s, int n, T const * vlist)
            {
                MatrixViewBase::Assign<T>(s, n, vlist);
            }

            //! @brief 深度拷贝 vlist 到 ridx, cidx 所指的子阵中
            //!
            //! @param [in] r 子阵起始行索引
            //! @param [in] c 子阵起始列索引
            //! @param [in] m 子阵行数
            //! @param [in] n 子阵列数
            //! @param [in] vlist 用户维护的内存
            void Assign(int r, int c, int m, int n, T const * vlist)
            {
                MatrixViewBase::Assign<T>(r, c, m, n, vlist);
            }

            //! @brief 交换 i, j 两行
            inline void RowSwap(int i, int j)
            {
                MatrixViewBase::RowSwap<T>(i, j);
            }

            //! @brief 交换 i, j 两列
            inline void ColSwap(int i, int j)
            {
                MatrixViewBase::ColSwap<T>(i, j);
            }


            friend std::ostream & operator << (std::ostream & s, MatrixView const & m)
            {
                s << std::endl;
                for (int ridx = 0; ridx < m.mNumRows; ridx++) {
                    s << m(ridx, 0);
                    for (int cidx = 1; cidx < m.mNumCols; cidx++)
                        s << "\t" << m(ridx, cidx);
                    s << std::endl;
                }
                return s;
            }

    };
}
}


#endif
