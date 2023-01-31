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

    template <typename T>
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
            //! @param [in] option 行优先?列优先?存储
            Matrix(int rows, int cols, EStorageOptions option = EStorageOptions::eColMajor)
            {
                __AllocDataBuffer__(rows * cols);
                __AllocView__(rows * cols);
                mView->Reshape(rows, cols, option);
            }


            //! @brief 构造函数
            //!
            //! @param [in] rows 矩阵行数
            //! @param [in] cols 矩阵列数
            //! @param [in] v 元素初值
            //! @param [in] option 行优先?列优先?存储
            Matrix(int rows, int cols, DataType const & v, EStorageOptions option = EStorageOptions::eColMajor)
            {
                Resize(rows, cols, option);
                Full(v);
            }

            //! @brief 析构函数
            ~Matrix()
            {
                if (EAlloc_User & mAlloclFlags)
                    return;

                __FreeView__();
                __FreeDataBuffer__();
            }

            //! @brief 改变矩阵尺寸，如果内存不够则重新分配。
            //!
            //! @note 如果矩阵数据是浅拷贝的则不能使用本函数
            //!
            //! @param [in] rows 新尺寸的行数
            //! @param [in] cols 新尺寸的列数
            //! @param [in] option 行优先?列优先?存储
            void Resize(int rows, int cols, EStorageOptions option = EStorageOptions::eColMajor)
            {
                int num = rows * cols;

                if (num <= 0) {
                    __FreeView__();
                    __FreeDataBuffer__();
                    return;
                }

                int storSize = mStorEnd - mStorBegin;
                if (storSize < num) {
                    if ((nullptr != mStorBegin) && !(EAlloc_DataBuffer & mAlloclFlags))
                        throw NotDeepCopyException("对象不是深度拷贝，需要 Clone()");
                    __AllocDataBuffer__(rows * cols);
                    __AllocView__(rows * cols);
                }

                mView->Reshape(rows, cols, option);
            }


            //! @brief 填充矩阵
            //!
            //! @param [in] v 目标值
            void Full(DataType const & v)
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
                    mStorBegin = nullptr;
                    mStorEnd = nullptr;

                    mAlloclFlags &= ~EAlloc_DataBuffer;
                }
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

            void __FreeView__()
            {
                if ((EAlloc_View & mAlloclFlags) && (nullptr != mView)) {
                    delete mView;
                    mView = nullptr;
                    mAlloclFlags &= ~EAlloc_View;
                }
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

            friend std::ostream & operator << (std::ostream & s, Matrix const & m)
            {
                return (nullptr == m.mView) ? s : (s << *(m.mView));
            }
    };

}
}


#endif
