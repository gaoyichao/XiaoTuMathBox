#ifndef XTMB_LA_MATRIX_DENSE_DMATRIX_H
#define XTMB_LA_MATRIX_DENSE_DMATRIX_H


#include <vector>
#include <iostream>
#include <initializer_list>


namespace xiaotu::math {

    template <typename _Scalar, EAlignType _align>
    struct Traits<DMatrix<_Scalar, _align>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreDyna;
    };

    /**
     * @brief 稠密矩阵
     * 
     * 以 std::vector<Scalar> 保存数据，适用与矩阵尺寸较大的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用小。
     * 矩阵尺寸小时，大量使用，将导致内存碎片化。
     */
    template <typename Scalar, EAlignType Align>
    class DMatrix : public MatrixBase<DMatrix<Scalar, Align>>
    {
        public:
            typedef MatrixBase<DMatrix> Base;

            using Base::NumDatas;
            using Base::At;
            using Base::Dot;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

            typedef DMatrixView<Scalar, Align> MatView;
            typedef DMatrixView<const Scalar, Align> CMatView;

        public:
            DMatrix(int rows, int cols)
                : mData(rows * cols), mRows(rows), mCols(cols)
            {}

            //! @brief 拷贝构造
            DMatrix(DMatrix const & mv)
                : mData(mv.mData), mRows(mv.Rows()), mCols(mv.Cols())
            {}

            //! @brief 拷贝赋值
            DMatrix & operator = (DMatrix const & mv)
            {
                mData = mv.mData;
                mRows = mv.mRows;
                mCols = mv.mCols;
                return *this;
            }

            //! @brief 构造一个全零矩阵
            static DMatrix Zero(int rows, int cols)
            {
                DMatrix re(rows, cols);
                re.Zeroing();
                return re;
            }

            //! @brief 构造一个单位矩阵
            static DMatrix Eye(int rows, int cols)
            {
                DMatrix re(rows, cols);
                re.Identity();
                return re;
            }

        public:
            //! @brief 转置
            DMatrix Transpose() const
            {
                DMatrix re(Cols(), Rows());
                xiaotu::math::Transpose(*this, re);
                return re;
            }

            //! @brief 获取视图
            inline MatView View() { return MatView(mData.data(), mRows, mCols); }
            //! @brief 获取视图
            inline CMatView View() const { return CMatView(mData.data(), mRows, mCols); }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mData.data(); }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mData.data(); }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mRows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mCols; }

        private:
            //! @brief 矩阵数据缓存
            std::vector<Scalar> mData;
            int mRows;
            int mCols;
    };


}

#endif


