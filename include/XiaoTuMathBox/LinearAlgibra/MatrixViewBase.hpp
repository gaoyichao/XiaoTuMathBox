#ifndef XTMB_LA_MATRIX_VIEW_BASE_H
#define XTMB_LA_MATRIX_VIEW_BASE_H

#include <stdint.h>
#include <sys/types.h>
#include <cassert>
#include <algorithm>

namespace xiaotu {
namespace math {

    //! @brief 矩阵视图的基类
    class MatrixViewBase {
        public:
            //! @brief 存储方式
            enum EStorageOptions {
                //! 列优先存储
                eColMajor = 0,
                //! 行优先存储
                eRowMajor = 1
            };
        public:

            //! @brief 构造函数
            //!
            //! @param [in] buffer 数据缓存
            //! @param [in] bytes 缓存字节数
            MatrixViewBase(uint8_t * buffer, size_t bytes)
            {
                mStorBegin = buffer;
                mBytes = bytes;

                mNumRows = 0;
                mNumCols = 0;
            }

            //! @brief 浅拷贝构造
            MatrixViewBase(MatrixViewBase const & mv)
            {
                mStorBegin = mv.mStorBegin;
                mBytes = mv.mBytes;

                Reshape(mv.mNumRows, mv.mNumCols, mv.mDataBytes, mv.mOption);
            }

            //! @brief 浅拷贝赋值
            MatrixViewBase & operator = (MatrixViewBase const & mv)
            {
                mStorBegin = mv.mStorBegin;
                mBytes = mv.mBytes;
                Reshape(mv.mNumRows, mv.mNumCols, mv.mDataBytes, mv.mOption);
                return *this;
            }

            //! @brief 重置矩阵形状，不重新分配内存
            //!
            //! @param [in] rows 新形状的行数
            //! @param [in] cols 新形状的列数
            //! @param [in] bytes 矩阵中每个元素的字节数 sizeof(T)
            //! @param [in] option 缓存的解析方式，默认列优先解析
            void Reshape(int rows, int cols, int bytes, EStorageOptions option = eColMajor)
            {
                assert(mBytes >= rows * cols * bytes);

                mNumRows = rows;
                mNumCols = cols;
                mDataBytes = bytes;
                mOption = option;
            }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                return (eColMajor == mOption) ? mNumRows * col + row : mNumCols * row + col;
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素地址
            template <typename T>
            inline T * Ptr(int idx)
            {
                int offset = mDataBytes * idx;
                assert(offset < mBytes);
                return (T*)(mStorBegin + offset);
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素地址
            template <typename T>
            inline T const * Ptr(int idx) const
            {
                int offset = mDataBytes * idx;
                assert(offset < mBytes);
                return (T*)(mStorBegin + offset);
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素地址
            template <typename T>
            inline T * Ptr(int row, int col)
            {
                assert(sizeof(T) == mDataBytes);
                return this->Ptr<T>(Idx(row, col));
            }

            //! @brief 获取指定位置的元素地址
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素地址
            template <typename T>
            inline T const * Ptr(int row, int col) const
            {
                assert(sizeof(T) == mDataBytes);
                return this->Ptr<T>(Idx(row, col));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            template <typename T>
            inline T & At(int idx)
            {
                return *(Ptr<T>(idx));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            template <typename T>
            inline T const & At(int idx) const
            {
                return *(Ptr<T>(idx));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            template <typename T>
            inline T & At(int row, int col)
            {
                return *(Ptr<T>(row, col));
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            template <typename T>
            inline T const & At(int row, int col) const
            {
                return *(Ptr<T>(row, col));
            }

            //! @brief 填充整个矩阵
            //!
            //! @param [in] 填充数值
            template <typename T>
            void Full(T const & v)
            {
                assert(sizeof(T) == mDataBytes);
                for (int ridx = 0; ridx < mNumRows; ridx++)
                    for (int cidx = 0; cidx < mNumCols; cidx++)
                        this->At<T>(ridx, cidx) = v;
            }

            //! @brief 深度拷贝 vlist 到 s 所指的连续内存中
            //!
            //! @param [in] s 目标起始索引
            //! @param [in] n 拷贝数据长度
            //! @param [in] vlist 用户维护的内存
            template <typename T>
            void Assign(int s, int n, T const * vlist)
            {
                assert((s + n) <= NumDatas());
                T * start = Ptr<T>(s);
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
            template <typename T>
            void Assign(int r, int c, int m, int n, T const * vlist)
            {
                int rend = r + m;
                int cend = c + n;

                assert(r >= 0 && rend <= mNumRows);
                assert(c >= 0 && cend <= mNumCols);

                int vidx = 0;
                for (int ridx = r; ridx < rend; ridx++)
                    for (int cidx = c; cidx < cend; cidx++)
                        At<T>(ridx, cidx) = vlist[vidx++];
            }

            //! @brief 交换 i, j 两行
            template <typename T>
            inline void RowSwap(int i, int j)
            {
                assert(i != j);
                assert(i < mNumRows && j < mNumRows);
                for (int k = 0; k < mNumCols; k++)
                    std::swap(At<T>(i, k), At<T>(j, k));
            }

            //! @brief 交换 i, j 两列
            template <typename T>
            inline void ColSwap(int i, int j)
            {
                assert(i != j);
                assert(i < mNumCols && j < mNumCols);
                for (int k = 0; k < mNumRows; k++)
                    std::swap(At<T>(k, i), At<T>(k, j));
            }
 
        public:
            //! @brief 获取缓存的字节数量
            size_t Bytes() const { return mBytes; }
            //! @brief 获取矩阵行数
            int Rows() const { return mNumRows; }
            //! @brief 获取矩阵列数
            int Cols() const { return mNumCols; }
            //! @brief 获取矩阵元素数量
            int NumDatas() const { return mNumCols * mNumRows; }
            //! @brief 获取矩阵元素字节数
            int DataBytes() const { return mDataBytes; }
            //! @brief 获取矩阵的存储方式
            EStorageOptions GetStorageOption() const { return mOption; }
        protected:
            //! @brief 矩阵数据起始地址
            uint8_t * mStorBegin = nullptr;
            //! @brief 缓存的字节数
            size_t mBytes;

            //! @brief 矩阵行数
            int mNumRows;
            //! @brief 矩阵列数
            int mNumCols;
            //! @brief 矩阵元素的字节数
            int mDataBytes;
            //! @brief 行优先?列优先?
            EStorageOptions mOption;
    };
}
}

#endif
