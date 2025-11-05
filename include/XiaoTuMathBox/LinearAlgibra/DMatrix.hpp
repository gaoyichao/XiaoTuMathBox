#ifndef XTMB_LA_MATRIX_DENSE_DMATRIX_H
#define XTMB_LA_MATRIX_DENSE_DMATRIX_H


#include <vector>
#include <iostream>
#include <initializer_list>


namespace xiaotu {

    template <typename _Scalar, EAlignType _align>
    struct Traits<DMatrix<_Scalar, _align>> {
        typedef _Scalar Scalar;
        constexpr static EAlignType Align = _align;
        constexpr static EStoreType Store = EStoreType::eStoreDyna;
    };

    template <typename _Scalar, EAlignType _align>
    struct Traits<const DMatrix<_Scalar, _align>> {
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

            /**
             * @brief 默认构造函数
             */
            DMatrix()
                : mRows(0), mCols(0)
            {}

            /**
             * @brief 构造函数
             * 
             * @param [in] rows 行数
             * @param [in] cols 列数
             */
            DMatrix(int rows, int cols)
                : mData(rows * cols), mRows(rows), mCols(cols)
            {}

            /**
             * @brief 拷贝构造(深度)
             * 
             * @param [in] mv 拷贝对象
             */
            DMatrix(DMatrix const & mv)
                : mData(mv.mData), mRows(mv.Rows()), mCols(mv.Cols())
            {}

            /**
             * @brief 拷贝构造(深度)
             * 
             * @param [in] mv 拷贝对象
             */
            template <typename Mat, bool IsMatrix = Mat::IsMatrix>
            DMatrix(Mat const & mv)
            {
                mRows = mv.Rows();
                mCols = mv.Cols();
                mData.resize(mRows * mCols);

                Assign(mv);
            }

            /**
             * @brief 拷贝赋值(深度)
             * 
             * @param [in] mv 拷贝对象
             */
            DMatrix & operator = (DMatrix const & mv)
            {
                mData = mv.mData;
                mRows = mv.mRows;
                mCols = mv.mCols;
                return *this;
            }

            /**
             * @brief 构造一个全零矩阵
             * 
             * @param [in] rows 行数
             * @param [in] cols 列数
             */
            static DMatrix Zero(int rows, int cols)
            {
                DMatrix re(rows, cols);
                re.Zeroing();
                return re;
            }

            /**
             * @brief 构造一个单位矩阵
             * 
             * @param [in] rows 行数
             * @param [in] cols 列数
             */
            static DMatrix Eye(int rows, int cols)
            {
                DMatrix re(rows, cols);
                re.Identity();
                return re;
            }

            /**
             * @brief 重新分配内存,不保存原始数据 
             * 
             * @param [in] rows 行数
             * @param [in] cols 列数
             */
            DMatrix & Resize(int rows, int cols)
            {
                mRows = rows;
                mCols = cols;
                mData.resize(mRows * mCols);
                return *this;
            }

        public:
            //! @brief 转置
            DMatrix Transpose() const
            {
                DMatrix re(Cols(), Rows());
                xiaotu::Transpose(*this, re);
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
            int mRows;
            int mCols;
    };


}

#endif


