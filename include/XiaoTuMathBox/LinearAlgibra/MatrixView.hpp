#ifndef XTMB_LA_MATRIX_VIEW_H
#define XTMB_LA_MATRIX_VIEW_H

#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename _Scalar, int _numRows, int _numCols, EStorageOptions _option>
    struct Traits<MatrixView<_Scalar, _numRows, _numCols, _option>> {
        typedef _Scalar Scalar;
        constexpr static int NumRows = _numRows;
        constexpr static int NumCols = _numCols;
        constexpr static int NumElements = _numRows * _numCols;
        constexpr static EStorageOptions StorOption = _option;
    };

    //! @brief 矩阵视图
    template <typename Scalar, int NumRows, int NumCols, EStorageOptions _option>
    class MatrixView : public MatrixBase<MatrixView<Scalar, NumRows, NumCols, _option>>
    {
        public:
            typedef MatrixBase<MatrixView> Base;

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
            //! @param [in] buffer 数据缓存
            //! @param [in] num 数据长度
            MatrixView(Scalar * buffer)
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
        private:
            //! @brief 矩阵数据起始地址
            Scalar * mStorBegin = nullptr;
    };

    //! @brief 列向量视图
    template <typename _Scalar, int _numRows, EStorageOptions _option = EStorageOptions::eColMajor>
    using VectorView = MatrixView<_Scalar, _numRows, 1, _option>;

    //! @brief 行向量视图
    template <typename _Scalar, int _numCols, EStorageOptions _option = EStorageOptions::eColMajor>
    using RowVectorView = MatrixView<_Scalar, 1, _numCols, _option>;

}


#endif
