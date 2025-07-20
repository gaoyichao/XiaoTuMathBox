#ifndef XTMB_LA_MATRIX_COL_VIEW_H
#define XTMB_LA_MATRIX_COL_VIEW_H

#include <cassert>
#include <iostream>
#include <vector>
#include <initializer_list>


namespace xiaotu::math {

    template <typename Derived>
    struct Traits<MatrixColView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    template <typename Derived>
    struct Traits<const MatrixColView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };


    //! @brief 选中特定列的视图
    template <typename Derived>
    class MatrixColView : public MatrixBase<MatrixColView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixColView> Base;
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
            //! @param [in] cols 选中的列索引, 将按照索引顺序输出
            MatrixColView(Derived & m, std::vector<int> const & col_list)
                : mMatrix(m), mSelectedCols(col_list)
            {}

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mMatrix.StorBegin(); }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mMatrix.Rows(); }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mSelectedCols.size(); }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < Rows() && col < Cols());
                int c = mSelectedCols[col];
                return mMatrix.Idx(row, c);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived & mMatrix;
            std::vector<int> mSelectedCols;
    };


    template <typename Derived>
    struct Traits<MatrixConstColView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    //! @brief 只读子矩阵视图
    template <typename Derived>
    class MatrixConstColView : public MatrixBase<MatrixConstColView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixConstColView> Base;
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
            //! @param [in] cols 选中的列索引, 将按照索引顺序输出
            MatrixConstColView(Derived const & m, std::vector<int> const & col_list)
                : mMatrix(m), mSelectedCols(col_list)
            {}

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mMatrix.Rows(); }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mSelectedCols.size(); }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < Rows() && col < Cols());
                int c = mSelectedCols[col];
                return mMatrix.Idx(row, c);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived const & mMatrix;
            std::vector<int> mSelectedCols;
    };

}


#endif
