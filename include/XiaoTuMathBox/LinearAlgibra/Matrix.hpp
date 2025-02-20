#ifndef XTMB_LA_MATRIX_DENSE_H
#define XTMB_LA_MATRIX_DENSE_H

#include <vector>
#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename _Scalar, int _numRows, int _numCols, EStorageOptions _option>
    struct Traits<Matrix<_Scalar, _numRows, _numCols, _option>> {
        typedef _Scalar Scalar;
        constexpr static int NumRows = _numRows;
        constexpr static int NumCols = _numCols;
        constexpr static int NumElements = _numRows * _numCols;
        constexpr static EStorageOptions StorOption = _option;
    };

    //! @brief 稠密矩阵
    template <typename Scalar, int NumRows, int NumCols, EStorageOptions StorOption>
    class Matrix : public MatrixBase<Matrix<Scalar, NumRows, NumCols, StorOption>>
    {
        public:
            typedef MatrixBase<Matrix> Base;

            using Base::At;
            using Base::operator();
            using Base::Assign;

            typedef MatrixView<Scalar, NumRows, NumCols, StorOption> MatView;
            typedef MatrixView<const Scalar, NumRows, NumCols, StorOption> CMatView;

        public:
            Matrix()
                : mData(NumRows * NumCols)
            {}

            //! @brief 拷贝构造
            Matrix(Matrix const & mv)
                : mData(mv.mData)
            {}

            //! @brief 构造,初始化列表
            Matrix(std::initializer_list<Scalar> && li)
                : mData(NumRows * NumCols)
            {
                MatView view = View();
                view = std::move(li);
            }

            //! @brief 拷贝赋值
            Matrix & operator = (Matrix const & mv)
            {
                mData = mv.mData;
                return *this;
            }

            //! @brief 构造一个全零矩阵
            static Matrix Zero()
            {
                Matrix re;
                re.Zeroing();
                return re;
            }

            //! @brief 构造一个单位矩阵
            static Matrix Eye()
            {
                Matrix re;
                re.Identity();
                return re;
            }

        public:
            //! @brief 转置
            Matrix<Scalar, NumCols, NumRows, StorOption> Transpose() const
            {
                Matrix<Scalar, NumCols, NumRows, StorOption> re;
                xiaotu::math::Transpose(*this, re);
                return re;
            }

            //! @brief 获取视图
            inline MatView View() { return MatView(mData.data()); }
            //! @brief 获取视图
            inline CMatView View() const { return CMatView(mData.data()); }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mData.data(); }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mData.data(); }

        private:
            //! @brief 矩阵数据缓存
            std::vector<Scalar> mData;
    };

    //! @brief 稠密列向量
    template <typename _Scalar, int _numRows, EStorageOptions _option = EStorageOptions::eColMajor>
    using Vector = Matrix<_Scalar, _numRows, 1, _option>;

    //! @brief 稠密行向量
    template <typename _Scalar, int _numCols, EStorageOptions _option = EStorageOptions::eColMajor>
    using RowVector = Matrix<_Scalar, 1, _numCols, _option>;
}


#endif
