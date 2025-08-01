#ifndef XTMB_LA_MATRIX_DENSE_H
#define XTMB_LA_MATRIX_DENSE_H

#include <vector>
#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreVector>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreVector;
    };

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<const Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreVector>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreVector;
    };

    /**
     * @brief 稠密矩阵
     * 
     * 以 std::vector<Scalar> 保存数据，适用与矩阵尺寸较大的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用小。
     * 矩阵尺寸小时，大量使用，将导致内存碎片化。
     */
    template <typename Scalar, int _rows, int _cols, EAlignType Align>
    class Matrix<Scalar, _rows, _cols, Align, EStoreType::eStoreVector> :
    public MatrixBase<Matrix<Scalar, _rows, _cols, Align>>
    {
        public:
            typedef MatrixBase<Matrix> Base;

            using Base::NumDatas;
            using Base::At;
            using Base::Dot;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

            constexpr static int NumRows = _rows;
            constexpr static int NumCols = _cols;

            typedef MatrixView<Scalar, _rows, _cols, Align> MatView;
            typedef MatrixView<const Scalar, _rows, _cols, Align> CMatView;

        public:
            Matrix()
                : mData(_rows * _cols)
            {}

            //! @brief 拷贝构造
            Matrix(Matrix const & mv)
                : mData(mv.mData)
            {}

            template <typename M>
            Matrix(M const & mv)
                : mData(_rows * _cols)
            {
                assert(Rows() == mv.Rows());
                assert(Cols() == mv.Cols());
                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = mv(ridx, cidx);
            }

            //! @brief 构造,初始化列表
            Matrix(std::initializer_list<Scalar> && li)
                : mData(_rows * _cols)
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
            Matrix<Scalar, _cols, _rows, Align, eStoreVector> Transpose() const
            {
                Matrix<Scalar, _cols, _rows, Align, eStoreVector> re;
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
                return (EAlignType::eRowMajor == Align)
                      ? Cols() * row + col
                      : Rows() * col + row;
            }

        private:
            //! @brief 矩阵数据缓存
            std::vector<Scalar> mData;
    };

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreArray>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreArray;
    };

    template <typename _Scalar, int _numRows, int _numCols, EAlignType _align>
    struct Traits<const Matrix<_Scalar, _numRows, _numCols, _align, EStoreType::eStoreArray>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreArray;
    };

    /**
     * @brief 稠密矩阵
     * 
     * 以 Scalar[] 保存数据，适用与矩阵尺寸较小的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用大。
     */
    template <typename Scalar, int _rows, int _cols, EAlignType Align>
    class Matrix<Scalar, _rows, _cols, Align, EStoreType::eStoreArray> :
    public MatrixBase<Matrix<Scalar, _rows, _cols, Align, EStoreType::eStoreArray>>
    {
        public:
            typedef MatrixBase<Matrix> Base;

            using Base::NumDatas;
            using Base::At;
            using Base::Dot;
            using Base::Assign;
            using Base::operator();
            using Base::operator=;

            constexpr static int NumRows = _rows;
            constexpr static int NumCols = _cols;

            typedef MatrixView<Scalar, _rows, _cols, Align> MatView;
            typedef MatrixView<const Scalar, _rows, _cols, Align> CMatView;

        public:
            Matrix()
            {}

            //! @brief 拷贝构造
            Matrix(Matrix const & mv)
            {
                int num = NumDatas();
                for (int i = 0; i < num; i++)
                    mData[i] = mv.mData[i];
            }

            template <typename M>
            Matrix(M const & mv)
            {
                assert(Rows() == mv.Rows());
                assert(Cols() == mv.Cols());
                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = mv(ridx, cidx);
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
            Matrix<Scalar, _cols, _rows, Align, eStoreArray> Transpose() const
            {
                Matrix<Scalar, _cols, _rows, Align, eStoreArray> re;
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
                return (EAlignType::eRowMajor == Align)
                      ? Cols() * row + col
                      : Rows() * col + row;
            }

        private:
            //! @brief 矩阵数据缓存
            Scalar mData[_rows * _cols];
    };




}


#endif
