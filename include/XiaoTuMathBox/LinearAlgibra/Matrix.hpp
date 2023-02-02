#ifndef XTMB_LA_MATRIX_H
#define XTMB_LA_MATRIX_H

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <initializer_list>

#include <XiaoTuMathBox/Exception.hpp>
#include <XiaoTuMathBox/LinearAlgibra/MatrixView.hpp>


namespace xiaotu {
namespace math {

    template <typename T, MatrixViewBase::EStorageOptions storMajor = MatrixViewBase::EStorageOptions::eColMajor>
    class Matrix {
        public:
            typedef T  DataType;
            typedef T* DataTypePtr;
            typedef MatrixViewBase::EStorageOptions EStorageOptions;

        public:
            //! @brief 默认构造函数
            Matrix()
            {}

            //! @brief 构造函数
            //!
            //! @param [in] rows 矩阵行数
            //! @param [in] cols 矩阵列数
            Matrix(int rows, int cols)
            {
                __AllocDataBuffer__(rows * cols);
                __AllocView__(rows * cols);
                mView->Reshape(rows, cols, storMajor);
            }

            //! @brief 构造函数
            //!
            //! @param [in] rows 矩阵行数
            //! @param [in] cols 矩阵列数
            //! @param [in] v 元素初值
            Matrix(int rows, int cols, DataType const & v)
            {
                Resize(rows, cols);
                Full(v);
            }

            //! @brief 根据初始化列表构造函数, 生成 mx1 的矩阵
            //!
            //! @param [in] li 初始化列表
            Matrix(std::initializer_list<DataType> li)
            {
                Resize(li.size(), 1);
                mView->Assign(0, li.size(), li.begin());
            }

            //! @brief 根据2D初始化列表构造函数
            //!
            //! @param [in] li 2D初始化列表
            Matrix(std::initializer_list<std::initializer_list<DataType>> li)
            {
                int rows = li.size();
                int cols = 0;

                for (auto rit = li.begin(); rit != li.end(); rit++)
                    cols = (rit->size() > cols) ? rit->size() : cols;

                Resize(rows, cols);
                mView->Assign(li);
            }

            //! @brief 构造函数, 深度拷贝 vlist
            //!
            //! @param [in] rows 矩阵行数
            //! @param [in] cols 矩阵列数
            //! @param [in] vlist 用户维护的内存
            Matrix(int rows, int cols, DataType const * vlist)
            {
                Resize(rows, cols);
                mView->Assign(0, 0, rows, cols, vlist);
            }

            //! @brief 浅拷贝构造
            Matrix(Matrix & mat)
            {
                mStorBegin = mat.mStorBegin;
                mStorEnd = mat.mStorEnd;
                mView = mat.mView;
            }

            //! @brief 深度拷贝构造
            Matrix(Matrix const & mat)
            {
                Resize(mat.Rows(), mat.Cols());
                mView->Assign(0, mat.NumDatas(), mat.mView->Ptr());
            }

            //! @brief 浅拷贝赋值
            Matrix & operator = (Matrix & mat)
            {
                Clear();
                mStorBegin = mat.mStorBegin;
                mStorEnd = mat.mStorEnd;
                mView = mat.mView;
                return *this;
            }

            //! @brief 深度拷贝赋值
            Matrix & operator = (Matrix const & mat)
            {
                Clear();
                Resize(mat.Rows(), mat.Cols());
                mView->Assign(0, mat.NumDatas(), mat.mView->Ptr());
                return *this;
            }


            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            Matrix & operator = (std::initializer_list<DataType> li)
            {
                mView->Assign(0, li.size(), li.begin());
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            //!
            //! 调用该函数之前，需要保证行数和列数能够容纳下 li 的所有数据。
            //! 不处理 li 没有覆盖到的部分。
            //!
            //! @param [in] li 用于拷贝的初始化列表
            Matrix & operator = (std::initializer_list<std::initializer_list<DataType>> li)
            {
                mView->Assign(li);
                return *this;
            }

            //! @brief 析构函数
            ~Matrix()
            {
                if (EAlloc_User & mAlloclFlags)
                    return;
                Clear();
            }

            inline void Clear()
            {
                __FreeView__();
                __FreeDataBuffer__();
            }

            //! @brief 改变矩阵尺寸，如果内存不够则重新分配。
            //!
            //! @note 如果矩阵数据是浅拷贝的则不能使用本函数
            //!
            //! @param [in] rows 新尺寸的行数
            //! @param [in] cols 新尺寸的列数
            void Resize(int rows, int cols)
            {
                int num = rows * cols;

                int storSize = mStorEnd - mStorBegin;
                if (storSize < num) {
                    if ((nullptr != mStorBegin) && !(EAlloc_DataBuffer & mAlloclFlags))
                        throw NotDeepCopyException("对象不是深度拷贝，需要 Clone()");
                    __AllocDataBuffer__(rows * cols);
                    __AllocView__(rows * cols);
                }

                if ((nullptr == mView) || !(EAlloc_View & mAlloclFlags))
                    __AllocView__(rows * cols);

                mView->Reshape(rows, cols, storMajor);
            }

            //! @brief 改变矩阵尺寸，如果内存不够则重新分配。
            //!
            //! @note 如果矩阵数据是浅拷贝的则不能使用本函数
            //!
            //! @param [in] rows 新尺寸的行数
            //! @param [in] cols 新尺寸的列数
            //! @param [in] v 填充的元素值
            inline void Resize(int rows, int cols, DataType const & v)
            {
                Resize(rows, cols);
                Full(v);
            }

            //! @brief 填充矩阵
            //!
            //! @param [in] v 目标值
            inline void Full(DataType const & v)
            {
                if (nullptr != mView)
                    mView->Full(v);
            }

        private:
            //! @brief 申请数据缓存
            //!
            //! @param [in] num 数据缓存大小
            void __AllocDataBuffer__(int num)
            {
                __FreeDataBuffer__();

                mStorBegin = new DataType[num];
                mStorEnd = mStorBegin + num;
                mAlloclFlags |= EAlloc_DataBuffer;
            }

            //! @brief 释放数据缓存
            void __FreeDataBuffer__()
            {
                if ((EAlloc_DataBuffer & mAlloclFlags) && (nullptr != mStorBegin)) {
                    delete [] mStorBegin;
                    mAlloclFlags &= ~EAlloc_DataBuffer;
                }
                mStorBegin = nullptr;
                mStorEnd = nullptr;
            }

            //! @brief 申请视图对象
            //!
            //! @param [in] num 矩阵元素数量
            void __AllocView__(int num)
            {
                assert(nullptr != mStorBegin);
                assert(num <= Storage());

                __FreeView__();
                mView = new MatrixView<T>(mStorBegin, num);
                mAlloclFlags |= EAlloc_View;
            }

            //! @brief 释放视图对象
            void __FreeView__()
            {
                if ((EAlloc_View & mAlloclFlags) && (nullptr != mView)) {
                    delete mView;
                    mAlloclFlags &= ~EAlloc_View;
                }
                mView = nullptr;
            }
 
        public:
            enum EAllocFlags {
                EAlloc_No         = 0x0000,   //! 没有申请内存
                EAlloc_DataBuffer = 0x0001,   //! 申请了数据缓存
                EAlloc_PtrBuffer  = 0x0002,   //! 申请了索引缓存
                EAlloc_View       = 0x0004,   //! 申请了视图对象
                EAlloc_All        = EAlloc_PtrBuffer | EAlloc_DataBuffer,
                EAlloc_User       = 0x1000,   //! 慎用!!! 析构时不考虑内存, 完全交给用户释放
            };

            //! @brief 获取内存维护标识
            inline int GetAllocFlags() const { return mAlloclFlags; }
            //! @brief 获取视图对象
            inline MatrixView<T> * GetView() { return mView; }
            //! @brief 获取视图对象
            inline MatrixView<T> const * GetView() const { return mView; }
        private:
            //! @brief 内存维护标识
            int mAlloclFlags = EAlloc_No;
            //! @brief 矩阵视图对象
            MatrixView<T> * mView = nullptr;

        public:
            //! brief 判定矩阵是否没有数据
            inline bool Empty() const { return mStorEnd == mStorBegin; }
            //! @brief 数据缓存长度
            inline int Storage() const { return mStorEnd - mStorBegin; }
        private:
            //! @brief 数据的起始地址
            DataType * mStorBegin = nullptr;
            //! @brief 数据的结束地址
            DataType * mStorEnd = nullptr;

        public:
            //! @brief 获取矩阵行数
            int Rows() const
            {
                return (nullptr == mView) ? 0 : mView->Rows();
            }

            //! @brief 获取矩阵列数
            int Cols() const
            {
                return (nullptr == mView) ? 0 : mView->Cols();
            }

            //! @brief 获取矩阵元素数量
            int NumDatas() const
            {
                return (nullptr == mView) ? 0 : mView->NumDatas();
            }

            //! @brief 交换 i, j 两行
            inline void RowSwap(int i, int j)
            {
                mView->RowSwap(i, j);
            }

            //! @brief 交换 i, j 两列
            inline void ColSwap(int i, int j)
            {
                mView->ColSwap(i, j);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline T & operator() (int row, int col)
            {
                return mView->At(row, col);
            }


            friend std::ostream & operator << (std::ostream & s, Matrix const & m)
            {
                return (nullptr == m.mView) ? s : (s << *(m.mView));
            }


    };

}
}


#endif
