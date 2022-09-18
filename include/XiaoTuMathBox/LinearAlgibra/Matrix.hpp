#ifndef XTMB_LA_MATRIX_H
#define XTMB_LA_MATRIX_H

#include <stdint.h>
#include <cassert>
#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>

namespace xiaotu {
namespace math {

    template <typename T>
    class Matrix {
        public:
            typedef T  DataType;
            typedef T* DataTypePtr;

        public:
            Matrix()
            {}

            //! @brief 构造函数
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            Matrix(int m, int n)
            {
                Resize(m, n);
            }

            //! @brief 构造函数
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            //! @param [in] v 元素初值
            Matrix(int m, int n, DataType const & v)
            {
                Resize(m, n);
                Full(v);
            }

            //! @brief 构造函数, 深度拷贝 vlist
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            //! @param [in] vlist 用户维护的内存
            Matrix(int m, int n, DataType const * vlist)
            {
                Resize(m, n);
                Assign(0, 0, m, n, vlist);
            }


            //! @brief 构造函数, 浅拷贝 vlist
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            //! @param [in] vlist 用户维护的内存
            Matrix(int m, int n, DataType * vlist)
            {
                mStorBegin = vlist;
                mStorEnd = mStorBegin + m * n;
                mAlloclFlags = EAlloc_User;
                __ReShape__(mStorBegin);
            }

            //! @brief 拷贝构造，深度
            Matrix(Matrix const & mat)
            {
                Resize(mat.mNumRows, mat.mNumCols);
                Assign(0, 0, mat.mNumRows, mat.mNumCols, mat.mStorBegin);
            }

            //! @brief 拷贝赋值，深度
            Matrix & operator = (Matrix const & mat)
            {
                Resize(mat.mNumRows, mat.mNumCols);
                Assign(0, 0, mat.mNumRows, mat.mNumCols, mat.mStorBegin);
                return *this;
            }

            //! @brief 析构函数
            ~Matrix()
            {
                __FreeData__();
                __FreeRowPtr__();
            }
        public:
            //! @brief 重置矩阵形状，不重新分配内存
            //!
            //! @param [in] m 新形状的行数
            //! @param [in] n 新形状的列数
            //!
            //! @return true 成功修改形状
            //!         false 新形状尺寸超出内存边界
            bool Reshape(int m, int n)
            {
                int storSize = mStorEnd - mStorBegin;
                int num = m * n;
                if (0 < num && num <= storSize) {
                    __ReShape__(m, n, mStorBegin);
                    return true;
                }

                return false;
            }
            //! @brief 改变矩阵尺寸，如果内存不够则重新分配。
            //!
            //! @note 如果矩阵数据是浅拷贝的则不能使用本函数
            //!
            //! @param [in] m 新尺寸的行数
            //! @param [in] m 新尺寸的列数
            void Resize(int m, int n)
            {
                assert(!(EAlloc_User & mAlloclFlags));

                int num = m * n;
                if (num <= 0) {
                    __FreeData__();
                    __FreeRowPtr__();
                    return;
                }
                int storSize = mStorEnd - mStorBegin;
                if (storSize < num) {
                    __FreeData__();

                    mStorBegin = new DataType[num];
                    mStorEnd = mStorBegin + num;
                    mAlloclFlags |= EAlloc_Data;
                }
                __ReShape__(m, n, mStorBegin);
            }
            //! @brief 改变矩阵尺寸并填充，如果内存不够则重新分配
            //!
            //! @param [in] m 新尺寸的行数
            //! @param [in] m 新尺寸的列数
            //! @param [in] v 填充的元素值
            //!
            //! @return false 如果矩阵数据是浅拷贝的则不能使用本函数
            //!         true  成功调整
            bool Resize(int m, int n, DataType const & v)
            {
                if (!Resize(m, n))
                    return false;

                Full(v);
                return true;
            }

            //! @brief 填充矩阵
            //!
            //! @param [in] v 目标值
            void Full(DataType const & v)
            {
                for (int ridx = 0; ridx < mNumRows; ridx++)
                    for (int cidx = 0; cidx < mNumCols; cidx++)
                        (*this)[ridx][cidx] = v;
            }

            //! @brief 深度拷贝 vlist 到 ridx, cidx 所指的子阵中
            //!
            //! @param [in] r 子阵起始行索引
            //! @param [in] c 子阵起始列索引
            //! @param [in] m 子阵行数
            //! @param [in] n 子阵列数
            //! @param [in] vlist 用户维护的内存
            void Assign(int r, int c, int m, int n, DataType const * vlist)
            {
                int rend = r + m;
                int cend = c + n;

                assert(r > 0 && rend < mNumRows);
                assert(c > 0 && cend < mNumCols);

                int vidx = 0;
                for (int ridx = r; ridx < rend; ridx++)
                    for (int cidx = c; cidx < cend; cidx++)
                        (*this)[ridx][cidx] = vlist[vidx++];
            }

            //! @brief 交换 i, j 两行
            inline void RowSwap(int i, int j)
            {
                assert(i != j);
                assert(i < mNumRows && j < mNumRows);
                for (int k = 0; k < mNumCols; k++)
                    std::swap((*this)[i][k], (*this)[j][k]);
            }

            //! @brief 交换 i, j 两列
            inline void ColSwap(int i, int j)
            {
                assert(i != j);
                assert(i < mNumCols && j < mNumCols);
                for (int k = 0; k < mNumRows; k++)
                    std::swap((*this)[k][i], (*this)[k][j]);
            }

        public:
            //! @brief 获取第 ridx 行的起始地址
            inline DataType * operator[] (int ridx)
            {
                assert(ridx < mNumRows);
                return mRowPtrBegin[ridx];
            }

            //! @brief 获取第 ridx 行的起始地址
            inline DataType const * operator[] (int ridx) const
            {
                assert(ridx < mNumRows);
                return mRowPtrBegin[ridx];
            }

            friend std::ostream & operator << (std::ostream & s, Matrix const & m)
            {
                s << std::endl;
                for (int ridx = 0; ridx < m.mNumRows; ridx++) {
                    s << m[ridx][0];
                    for (int cidx = 1; cidx < m.mNumCols; cidx++)
                        s << "\t" << m[ridx][cidx];
                    s << std::endl;
                }

                return s;
            }

        public:
            inline bool Empty() const { return mStorEnd == mStorBegin; }
            inline int Storage() const { return mStorEnd - mStorBegin; }
            inline int Elements() const { return mNumRows * mNumCols; }
            inline int Rows() const { return mNumRows; }
            inline int Cols() const { return mNumCols; }

        private:
            
            void __ReShape__(int m, int n, DataType * vlist)
            {
                int rowPtrSize = mRowPtrEnd - mRowPtrBegin;
                if (rowPtrSize < m) {
                    __FreeRowPtr__();

                    mRowPtrBegin = new DataTypePtr[m];
                    mRowPtrEnd = mRowPtrBegin + m;
                    mAlloclFlags |= EAlloc_RowPtr;
                }

                mNumRows = m;
                mNumCols = n;

                mRowPtrBegin[0] = vlist;
                for (int i = 1; i < m; i++)
                    mRowPtrBegin[i] = mRowPtrBegin[i-1] + mNumCols;
            }

            void __FreeRowPtr__()
            {
                if (EAlloc_RowPtr & mAlloclFlags) {
                    delete [] mRowPtrBegin;
                    mRowPtrBegin = nullptr;
                    mRowPtrEnd = nullptr;

                    mAlloclFlags &= ~EAlloc_RowPtr;
                }
            }

            void __FreeData__()
            {
                if (EAlloc_Data & mAlloclFlags) {
                    delete [] mStorBegin;
                    mStorBegin = nullptr;
                    mStorEnd = nullptr;

                    mAlloclFlags &= ~EAlloc_Data;
                }
            }

        public:
            enum EAllocFlags {
                EAlloc_No     = 0x00,   //! 没有申请内存
                EAlloc_RowPtr = 0x01,   //! 申请了行地址索引
                EAlloc_Data   = 0x02,   //! 申请了数据地址
                EAlloc_All    = EAlloc_RowPtr | EAlloc_Data,
                EAlloc_User   = 0x04,   //! 不维护内存
            };

            int GetAllocFlags() const { return mAlloclFlags; }
        private:
            int mAlloclFlags = EAlloc_No;
            int mNumRows = 0;
            int mNumCols = 0;

            DataType * mStorBegin = nullptr;
            DataType * mStorEnd = nullptr;
            DataTypePtr * mRowPtrBegin = nullptr;
            DataTypePtr * mRowPtrEnd = nullptr;
    };

}
}


#endif
