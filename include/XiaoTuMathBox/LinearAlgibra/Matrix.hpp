#ifndef XTMB_LA_MATRIX_DENSE_H
#define XTMB_LA_MATRIX_DENSE_H

#include <vector>
#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreVector>> {
        typedef _Scalar Scalar;
        constexpr static int NumRows = _numRows;
        constexpr static int NumCols = _numCols;
        constexpr static int NumElements = _numRows * _numCols;
        constexpr static EAlignType Align = _align;
    };

    /**
     * @brief 稠密矩阵
     * 
     * 以 std::vector<Scalar> 保存数据，适用与矩阵尺寸较大的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用小。
     * 矩阵尺寸小时，大量使用，将导致内存碎片化。
     */
    template <typename Scalar, int NumRows, int NumCols, EAlignType Align>
    class Matrix<Scalar, NumRows, NumCols, Align, EStoreType::eStoreVector> :
    public MatrixBase<Matrix<Scalar, NumRows, NumCols, Align>>
    {
        public:
            typedef MatrixBase<Matrix> Base;

            using Base::At;
            using Base::operator();
            using Base::Assign;

            typedef MatrixView<Scalar, NumRows, NumCols, Align> MatView;
            typedef MatrixView<const Scalar, NumRows, NumCols, Align> CMatView;

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
            Matrix<Scalar, NumCols, NumRows, Align, eStoreVector> Transpose() const
            {
                Matrix<Scalar, NumCols, NumRows, Align, eStoreVector> re;
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

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreArray>> {
        typedef _Scalar Scalar;
        constexpr static int NumRows = _numRows;
        constexpr static int NumCols = _numCols;
        constexpr static int NumElements = _numRows * _numCols;
        constexpr static EAlignType Align = _align;
    };
    /**
     * @brief 稠密矩阵
     * 
     * 以 Scalar[] 保存数据，适用与矩阵尺寸较小的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用大。
     */
    template <typename Scalar, int NumRows, int NumCols, EAlignType Align>
    class Matrix<Scalar, NumRows, NumCols, Align, EStoreType::eStoreArray> :
    public MatrixBase<Matrix<Scalar, NumRows, NumCols, Align, EStoreType::eStoreArray>>
    {
        public:
            typedef MatrixBase<Matrix> Base;

            using Base::NumDatas;
            using Base::At;
            using Base::operator();
            using Base::Assign;

            typedef MatrixView<Scalar, NumRows, NumCols, Align> MatView;
            typedef MatrixView<const Scalar, NumRows, NumCols, Align> CMatView;

        public:
            Matrix()
            {}

            //! @brief 拷贝构造
            Matrix(Matrix const & mv)
                : mData(mv.mData)
            {
                int num = NumDatas();
                for (int i = 0; i < num; i++)
                    mData[i] = mv.mData[i];
            }

            //! @brief 构造,初始化列表
            Matrix(std::initializer_list<Scalar> && li)
            {
                MatView view = View();
                view = std::move(li);
            }

            //! @brief 拷贝赋值
            Matrix & operator = (Matrix const & mv)
            {
                int num = NumDatas();
                for (int i = 0; i < num; i++)
                    mData[i] = mv.mData[i];
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
            Matrix<Scalar, NumCols, NumRows, Align, eStoreArray> Transpose() const
            {
                Matrix<Scalar, NumCols, NumRows, Align, eStoreArray> re;
                xiaotu::math::Transpose(*this, re);
                return re;
            }

            //! @brief 获取视图
            inline MatView View() { return MatView(mData); }
            //! @brief 获取视图
            inline CMatView View() const { return CMatView(mData); }

        public:
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar * StorBegin() { return mData; }
            //! @brief 获取矩阵数据存储的起始地址
            inline Scalar const * StorBegin() const { return mData; }

        private:
            //! @brief 矩阵数据缓存
            Scalar mData[NumRows * NumCols];
    };




}


#endif
