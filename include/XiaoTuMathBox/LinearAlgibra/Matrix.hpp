#ifndef XTMB_LA_MATRIX_H
#define XTMB_LA_MATRIX_H

#include <stdint.h>
#include <cassert>

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
                __Resize__(m, n);
            }

            //! @brief 构造函数
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            //! @param [in] v 元素初值
            Matrix(int m, int n, DataType const & v)
            {
                __Resize__(m, n);
                Full(v);
            }

            //! @brief 构造函数, 深度拷贝 vlist
            //!
            //! @param [in] m 矩阵行数
            //! @param [in] n 矩阵列数
            //! @param [in] vlist 用户维护的内存
            Matrix(int m, int n, DataType const * vlist)
            {
                __Resize__(m, n);
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
                __Resize__(mat.mNumRows, mat.mNumCols);
                Assign(0, 0, mat.mNumRows, mat.mNumCols, mat.mStorBegin);
            }

            //! @brief 拷贝赋值，深度
            Matrix & operator = (Matrix const & mat)
            {
                __Resize__(mat.mNumRows, mat.mNumCols);
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
            //! @brief 获取第 ridx 行的起始地址
            inline DataType * operator[] (int ridx)
            {
                assert(ridx < mNumRows);
                if (Empty())
                    return nullptr;
                return mRowPtrBegin[ridx];
            }

            //! @brief 获取第 ridx 行的起始地址
            inline DataType const * operator[] (int ridx) const
            {
                assert(ridx < mNumRows);
                return mRowPtrBegin[ridx];
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

        public:
            inline bool Empty() const { return mStorEnd == mStorBegin; }
            inline int Rows() const { return mNumRows; }
            inline int Cols() const { return mNumCols; }

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

            bool Resize(int m, int n)
            {
                if (EAlloc_User == mAlloclFlags)
                    return false;
                __Resize__(m, n);
                return true;
            }

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

            void __Resize__(int m, int n)
            {
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

        private:
            enum EAllocFlags {
                EAlloc_No     = 0x00,   //! 没有申请内存
                EAlloc_RowPtr = 0x01,   //! 申请了行地址索引
                EAlloc_Data   = 0x02,   //! 申请了数据地址
                EAlloc_All    = EAlloc_RowPtr | EAlloc_Data,
                EAlloc_User   = 0x04,   //! 不维护内存
            };
            EAllocFlags mAlloclFlags = EAlloc_No;
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
