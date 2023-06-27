#ifndef XTMB_LA_MATRIX_VIEW_H
#define XTMB_LA_MATRIX_VIEW_H

#include <iostream>
#include <initializer_list>

namespace xiaotu {
namespace math {

    template<typename T, int _Rows, int _Cols, EStorageOptions _Options>
    struct traits<MatrixView<T, _Rows, _Cols, _Options> >
    {
        public:
          typedef T Scalar;
    };
    
    //! @brief 矩阵视图
    template <typename T, int numRows, int numCols, EStorageOptions _option>
    class MatrixView : public MatrixViewBase<MatrixView<T, numRows, numCols, _option> >
    {
        public:
            typedef MatrixViewBase<MatrixView> Base;

        public:
            //! @brief 构造函数
            //!
            //! @param [in] buffer 数据缓存
            //! @param [in] num 数据长度
            MatrixView(T * buffer)
                : mStorBegin(buffer)
            {}

            //! @brief 浅拷贝构造
            MatrixView(MatrixView const & mv)
                : mStorBegin(mv.mStorBegin)
            {}

            //! @brief 浅拷贝赋值
            MatrixView & operator = (MatrixView const & mv)
            {
                mStorBegin = mv.mStorBegin;
                return *this;
            }

            //! @brief 浅拷贝赋值
            MatrixView & operator = (T * buffer)
            {
                mStorBegin = buffer;
                return *this;
            }
        public:
            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            MatrixView & operator = (std::initializer_list<T> li)
            {
                assert(li.size() == Base::NumDatas());
                Base::Assign(0, 0, Rows(), Cols(), li.begin());
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            //!
            //! 调用该函数之前，需要保证行数和列数能够容纳下 li 的所有数据。
            //! 不处理 li 没有覆盖到的部分。
            //!
            //! @param [in] li 用于拷贝的初始化列表
            MatrixView & operator = (std::initializer_list<std::initializer_list<T>> li)
            {
                Base::Assign(li);
                return *this;
            }

        public:
            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline virtual int Idx(int row, int col) const override
            {
                return (EStorageOptions::eColMajor == _option) ? Rows() * col + row : Cols() * row + col;
            }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline virtual T * StorBegin() override { return mStorBegin; }
            //! @brief 获取矩阵数据存储的起始地址
            inline virtual T const * StorBegin() const override { return mStorBegin; }
            //! @brief 获取矩阵行数
            inline virtual int Rows() const override { return numRows; }
            //! @brief 获取矩阵列数
            inline virtual int Cols() const override { return numCols; }
            //! @brief 获取矩阵的存储方式
            inline EStorageOptions GetStorageOption() const { return _option; }
        private:
            //! @brief 矩阵数据起始地址
            T * mStorBegin = nullptr;
    };

    template <typename T, EStorageOptions _option>
    class MatrixView<T, Dynamic, Dynamic, _option> : public MatrixViewBase<MatrixView<T, Dynamic, Dynamic, _option> >
    {
        public:
            typedef MatrixViewBase<MatrixView> Base;
 
        public:
            MatrixView(T * buffer, int cols, int rows)
                : mStorBegin(buffer), mCols(cols), mRows(rows)
            { }

            //! @brief 浅拷贝构造
            MatrixView(MatrixView const & mv)
                : mStorBegin(mv.mStorBegin), mCols(mv.mCols), mRows(mv.mRows)
            {}

            //! @brief 浅拷贝赋值
            MatrixView & operator = (MatrixView const & mv)
            {
                mStorBegin = mv.mStorBegin;
                mCols = mv.mCols;
                mRows = mv.Rows;
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            MatrixView & operator = (std::initializer_list<T> li)
            {
                assert(li.size() == Base::NumDatas());
                Base::Assign(0, 0, Rows(), Cols(), li.begin());
                return *this;
            }

            //! @brief 拷贝赋值，不改变矩阵形状，需要保证内存足够
            //!
            //! 调用该函数之前，需要保证行数和列数能够容纳下 li 的所有数据。
            //! 不处理 li 没有覆盖到的部分。
            //!
            //! @param [in] li 用于拷贝的初始化列表
            MatrixView & operator = (std::initializer_list<std::initializer_list<T>> li)
            {
                Base::Assign(li);
                return *this;
            }

            //! @brief 重置矩阵形状，不重新分配内存
            //!
            //! @param [in] rows 新形状的行数
            //! @param [in] cols 新形状的列数
            void Reshape(int rows, int cols)
            {
                assert(Base::NumDatas() >= rows * cols);

                mRows = rows;
                mCols = cols;
            }
        public:
            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline virtual int Idx(int row, int col) const override
            {
                return (EStorageOptions::eColMajor == _option) ? Rows() * col + row : Cols() * row + col;
            }


        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline virtual T * StorBegin() override { return mStorBegin; }
            //! @brief 获取矩阵数据存储的起始地址
            inline virtual T const * StorBegin() const override { return mStorBegin; }
            //! @brief 获取矩阵行数
            inline virtual int Rows() const override { return mRows; }
            //! @brief 获取矩阵列数
            inline virtual int Cols() const override { return mCols; }
            //! @brief 获取矩阵的存储方式
            inline EStorageOptions GetStorageOption() const { return _option; }
        private:
            //! @brief 矩阵数据起始地址
            T * mStorBegin = nullptr;
            int mRows;
            int mCols;
    };
}
}


#endif
