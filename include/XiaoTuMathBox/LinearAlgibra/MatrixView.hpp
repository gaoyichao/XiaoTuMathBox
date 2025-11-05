#ifndef XTMB_LA_MATRIX_VIEW_H
#define XTMB_LA_MATRIX_VIEW_H

#include <cassert>
#include <iostream>
#include <initializer_list>

namespace xiaotu {

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<MatrixView<_Scalar, _numRows, _numCols, _align>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreNone;
    };
    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<const MatrixView<_Scalar, _numRows, _numCols, _align>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreNone;
    };


    //! @brief 矩阵视图
    template <typename Scalar, int _rows, int _cols, EAlignType _align>
    class MatrixView : public MatrixBase<MatrixView<Scalar, _rows, _cols, _align>>
    {
        public:
            typedef MatrixBase<MatrixView> Base;

            using Base::At;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

            constexpr static int NumRows = _rows;
            constexpr static int NumCols = _cols;

        public:
            ////////////////////////////////////////////////////////
            //
            //  构造函数和拷贝
            //
            ////////////////////////////////////////////////////////

            //! @brief 构造函数
            //!
            //! @param [in] buffer 数据缓存
            //! @param [in] num 数据长度
            MatrixView(Scalar * buffer)
                : mStorBegin(buffer)
            {}

            //! @brief 浅拷贝构造
            MatrixView(MatrixView const & mv)
                : mStorBegin(mv.mStorBegin)
            {}

            //! @brief 浅拷贝构造
            template <typename M>
            MatrixView(M & mv)
                : mStorBegin(mv.StorBegin())
            {
                assert(_rows == mv.Rows());
                assert(_cols == mv.Cols());
            }

            //! @brief 浅拷贝赋值
            template <typename M>
            MatrixView & operator = (M const & mv)
            {
                mStorBegin = mv.StorBegin();
                return *this;
            }

            //! @brief 浅拷贝赋值
            MatrixView & operator = (Scalar * buffer)
            {
                mStorBegin = buffer;
                return *this;
            }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mStorBegin; }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mStorBegin; }

            //! @brief 获取矩阵行数
            inline int Rows() const { return _rows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return _cols; }

            //! @brief 计算指定行列索引的展开索引
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素的展开索引
            inline int Idx(int row, int col) const
            {
                return (EAlignType::eRowMajor == _align)
                      ? Cols() * row + col
                      : Rows() * col + row;
            }

        private:
            //! @brief 矩阵数据起始地址
            Scalar * mStorBegin = nullptr;
    };


}


#endif
