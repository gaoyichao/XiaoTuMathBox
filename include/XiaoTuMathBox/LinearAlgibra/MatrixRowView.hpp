#ifndef XTMB_LA_MATRIX_ROW_VIEW_H
#define XTMB_LA_MATRIX_ROW_VIEW_H

#include <cassert>
#include <iostream>
#include <vector>
#include <initializer_list>

namespace xiaotu {

    template <typename Derived>
    struct Traits<MatrixRowView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    template <typename Derived>
    struct Traits<const MatrixRowView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    //! @brief 选中特定行的视图
    template <typename Derived>
    class MatrixRowView : public MatrixBase<MatrixRowView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixRowView> Base;
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
            //! @param [in] cols 选中的行索引, 将按照索引顺序输出
            MatrixRowView(Derived & m, std::vector<int> const & row_list)
                : mMatrix(m), mSelectedRows(row_list)
            {}

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mMatrix.StorBegin(); }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mSelectedRows.size(); }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mMatrix.Cols(); }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < Rows() && col < Cols());
                int r = mSelectedRows[row];
                return mMatrix.Idx(r, col);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived & mMatrix;
            std::vector<int> mSelectedRows;
    };

    template <typename Derived>
    struct Traits<MatrixConstRowView<Derived>> {
        typedef typename Traits<Derived>::Scalar Scalar;
        constexpr static EAlignType Align = Derived::Align;
        constexpr static EStoreType Store = EStoreType::eStoreProxy;
    };

    //! @brief 只读子矩阵视图
    template <typename Derived>
    class MatrixConstRowView : public MatrixBase<MatrixConstRowView<Derived>>
    {
        public:
            typedef MatrixBase<MatrixConstRowView> Base;
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
            //! @param [in] cols 选中的行索引, 将按照索引顺序输出
            MatrixConstRowView(Derived const & m, std::vector<int> const & row_list)
                : mMatrix(m), mSelectedRows(row_list)
            {}

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mMatrix.StorBegin(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mSelectedRows.size(); }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mMatrix.Cols(); }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                assert(row < Rows() && col < Cols());
                int r = mSelectedRows[row];
                return mMatrix.Idx(r, col);
            }

            int MatrixRows() const { return mMatrix.Rows(); }
            int MatrixCols() const { return mMatrix.Cols(); }

        private:
            Derived const & mMatrix;
            std::vector<int> mSelectedRows;
    };


}

#endif
