#ifndef XTMB_LA_MATRIX_VIEW_D_H
#define XTMB_LA_MATRIX_VIEW_D_H

#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename T, EAlignType _align>
    struct Traits<DMatrixView<T, _align>> {
        typedef T Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreNone;
    };

    //! @brief 矩阵视图
    template <typename Scalar, EAlignType _align>
    class DMatrixView : public MatrixBase<DMatrixView<Scalar, _align>>
    {
        public:
            typedef MatrixBase<DMatrixView> Base;

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
            
            DMatrixView(Scalar * buffer, int rows, int cols)
                : mStorBegin(buffer), mRows(rows), mCols(cols)
            {}

            //! @brief 浅拷贝构造
            DMatrixView(DMatrixView const & mv)
                : mStorBegin(mv.mStorBegin), mRows(mv.Rows()), mCols(mv.Cols())
            {}

            //! @brief 浅拷贝构造
            template <typename M>
            DMatrixView(M & mv)
                : mStorBegin(mv.StorBegin()), mRows(mv.Rows()), mCols(mv.Cols())
            {}

            //! @brief 浅拷贝赋值
            template <typename M>
            DMatrixView & operator = (M & mv)
            {
                mStorBegin = mv.StorBegin();
                mRows = mv.Rows();
                mCols = mv.Cols();
                return *this;
            }

            //! @brief 浅拷贝赋值
            //! 只改变数据内容，不改变矩阵尺寸
            DMatrixView & operator = (Scalar * buffer)
            {
                mStorBegin = buffer;
                return *this;
            }

            //! @brief 修改矩阵尺寸
            //! @param [in] nRows 新矩阵行数
            //! @param [in] nCols 新矩阵列数
            void Reshape(int nRows, int nCols)
            {
                mRows = nRows;
                mCols = nCols;
            }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mStorBegin; }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mStorBegin; }

            //! @brief 获取矩阵行数
            inline int Rows() const { return mRows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return mCols; }
        private:
            //! @brief 矩阵数据起始地址
            Scalar * mStorBegin = nullptr;
            int mRows = 0;
            int mCols = 0;
    };


}


#endif
