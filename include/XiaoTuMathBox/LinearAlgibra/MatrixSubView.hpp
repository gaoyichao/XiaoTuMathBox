#ifndef XTMB_LA_MATRIX_SUB_VIEW_H
#define XTMB_LA_MATRIX_SUB_VIEW_H

#include <cassert>
#include <iostream>
#include <initializer_list>


namespace xiaotu::math {

    template <typename Derived>
    struct Traits<MatrixSubView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    //! @brief 子矩阵视图
    template <typename Derived>
    class MatrixSubView : public MatrixBase<MatrixSubView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixSubView> Base;
            typedef typename Traits<Derived>::Scalar Scalar;

            using Base::At;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

        public:
            ////////////////////////////////////////////////////////
            //
            //  构造函数和拷贝
            //
            ////////////////////////////////////////////////////////


            //! @brief 构造函数
            //!
            //! @param [in] m 目标矩阵
            //! @param [in] r 子阵的起始行
            //! @param [in] c 子阵的起始列
            //! @param [in] rows 子阵的行数
            //! @param [in] cols 子阵的列数
            MatrixSubView(Derived & m, int r, int c, int rows, int cols)
                : mMatrix(m), mStartRow(r), mStartCol(c), mRows(rows), mCols(cols)
            {}

            //! @brief 获取子阵, 相对于自己
            //!
            //! @param [in] r 子阵的起始行
            //! @param [in] c 子阵的起始列
            //! @param [in] rows 子阵的行数
            //! @param [in] cols 子阵的列数
            MatrixSubView SubMatrix(int r, int c, int rows, int cols)
            {
                return MatrixSubView(mMatrix, r + mStartRow, c + mStartCol, rows, cols);
            }

            //! @brief 获取行向量
            //!
            //! @param [in] r 目标行
            MatrixSubView Row(int r)
            {
                return MatrixSubView(mMatrix, r + mStartRow, mStartCol, 1, Cols());
            }

            //! @brief 获取列向量
            //!
            //! @param [in] c 目标行
            MatrixSubView Col(int c)
            {
                return MatrixSubView(mMatrix, mStartRow, c + mStartCol, Rows(), 1);
            }


            //! @brief 获取子阵
            //!
            //! @param [in] r 子阵的起始行
            //! @param [in] c 子阵的起始列
            //! @param [in] rows 子阵的行数
            //! @param [in] cols 子阵的列数
            MatrixConstSubView<Derived>
            SubMatrix(int r, int c, int rows, int cols) const
            {
                return MatrixConstSubView<Derived>(mMatrix, r + mStartRow, c + mStartCol, rows, cols);
            }

            //! @brief 获取行向量
            //!
            //! @param [in] r 目标行
            MatrixConstSubView<Derived>
            Row(int r) const
            {
                return MatrixConstSubView<Derived>(mMatrix, r + mStartRow, mStartCol, 1, Cols());
            }

            //! @brief 获取列向量
            //!
            //! @param [in] c 目标行
            MatrixConstSubView<Derived>
            Col(int c) const
            {
                return MatrixConstSubView<Derived>(mMatrix, mStartRow, c + mStartCol, Rows(), 1);
            }


        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mMatrix.StorBegin(); }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mRows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mCols; }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < mRows && col < mCols);
                int r = mStartRow + row;
                int c = mStartCol + col;
                return mMatrix.Idx(r, c);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived & mMatrix;
            int mStartRow;
            int mStartCol;
            int mRows;
            int mCols;
    };

    template <typename Derived>
    struct Traits<MatrixConstSubView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    //! @brief 只读子矩阵视图
    template <typename Derived>
    class MatrixConstSubView : public MatrixBase<MatrixConstSubView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixConstSubView> Base;
            typedef typename Traits<Derived>::Scalar Scalar;

            using Base::At;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

        public:
            ////////////////////////////////////////////////////////
            //
            //  构造函数和拷贝
            //
            ////////////////////////////////////////////////////////


            //! @brief 构造函数
            //!
            //! @param [in] m 目标矩阵
            //! @param [in] r 子阵的起始行
            //! @param [in] c 子阵的起始列
            //! @param [in] rows 子阵的行数
            //! @param [in] cols 子阵的列数
            MatrixConstSubView(Derived const & m, int r, int c, int rows, int cols)
                : mMatrix(m), mStartRow(r), mStartCol(c), mRows(rows), mCols(cols)
            {}

            //! @brief 获取子阵
            //!
            //! @param [in] r 子阵的起始行
            //! @param [in] c 子阵的起始列
            //! @param [in] rows 子阵的行数
            //! @param [in] cols 子阵的列数
            MatrixConstSubView SubMatrix(int r, int c, int rows, int cols) const
            {
                return MatrixConstSubView(mMatrix, r + mStartRow, c + mStartCol, rows, cols);
            }

            //! @brief 获取行向量
            //!
            //! @param [in] r 目标行
            MatrixConstSubView Row(int r) const
            {
                return MatrixConstSubView(mMatrix, r + mStartRow, mStartCol, 1, Cols());
            }

            //! @brief 获取列向量
            //!
            //! @param [in] c 目标行
            MatrixConstSubView Col(int c) const
            {
                return MatrixConstSubView(mMatrix, mStartRow, c + mStartCol, Rows(), 1);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar const & operator() (int idx)
            {
                return __At__(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar const & operator() (int row, int col)
            {
                return __At__(row, col);
            }

        private:

            inline Scalar const & __At__(int idx) const
            {
                return At(idx);
            }

            inline Scalar const & __At__(int row, int col) const
            {
                return At(row, col);
            }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mRows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mCols; }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < mRows && col < mCols);
                int r = mStartRow + row;
                int c = mStartCol + col;
                return mMatrix.Idx(r, c);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived const & mMatrix;
            int mStartRow;
            int mStartCol;
            int mRows;
            int mCols;
    };
}


#endif

