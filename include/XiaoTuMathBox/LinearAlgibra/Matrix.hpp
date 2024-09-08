#ifndef XTMB_LA_MATRIX_DENSE_H
#define XTMB_LA_MATRIX_DENSE_H

#include <vector>
#include <iostream>
#include <initializer_list>

namespace xiaotu::math {

     //! @brief 稠密矩阵
    template <typename _Scalar, int _numRows, int _numCols, EStorageOptions _option>
    class Matrix
    {
        public:
            typedef _Scalar Scalar;
            constexpr static int NumRows = _numRows;
            constexpr static int NumCols = _numCols;
            constexpr static EStorageOptions StorOption = _option;

            typedef MatrixView<Scalar, NumRows, NumCols, StorOption> MatView;
        public:
            Matrix()
                : mData(_numRows * _numCols), mView(mData.data())
            {}

            //! @brief 拷贝构造
            Matrix(Matrix const & mv)
                : mData(mv.mData), mView(mData.data())
            {}

            //! @brief 拷贝赋值
            Matrix & operator = (Matrix const & mv)
            {
                mData = mv.mData;
                mView = mData.data();
                return *this;
            }

            //! @brief 构造,初始化列表
            Matrix(std::initializer_list<Scalar> li)
                : mData(_numRows * _numCols), mView(mData.data())
            {
                mView = li;
            }

            //! @brief 深度拷贝 m
            //!
            //! @param [in] m 将矩阵 m 中的值拷贝过来
            void Assign(Matrix const & m)
            {
                assert(Rows() == m.Rows() && Cols() == m.Cols());

                for (int ridx = 0; ridx < Rows(); ridx++)
                    for (int cidx = 0; cidx < Cols(); cidx++)
                        At(ridx, cidx) = m(ridx, cidx);
            }

        public:
            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar & At(int idx)
            {
                return mView.At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar const & At(int idx) const
            {
                return mView.At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar & At(int row, int col)
            {
                return mView.At(row, col);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar const & At(int row, int col) const
            {
                return mView.At(row, col);
            }


            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar & operator() (int idx)
            {
                return mView.At(idx);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] idx 元素的展开索引
            //! @return 元素引用
            inline Scalar const & operator() (int idx) const
            {
                return mView.At(idx);
            }


            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar & operator() (int row, int col)
            {
                return mView.At(row, col);
            }

            //! @brief 获取指定位置的元素引用
            //!
            //! @param [in] row 行索引
            //! @param [in] col 列索引
            //! @return 元素引用
            inline Scalar const & operator() (int row, int col) const
            {
                return mView.At(row, col);
            }

        public:

            //! @brief 单位化, 将矩阵改写成单位矩阵
            void Identity()
            {
                mView.Identity();
            }

            //! @brief 交换 i, j 两行
            inline void RowSwap(int i, int j)
            {
                mView.RowSwap(i, j);
            }

            //! @brief 交换 i, j 两列
            inline void ColSwap(int i, int j)
            {
                mView.ColSwap(i, j);
            }

            inline MatView & View() { return mView; }
            inline MatView const & View() const { return mView; }

            friend std::ostream & operator << (std::ostream & s, Matrix const & m)
            {
                s << std::endl;
                for (int ridx = 0; ridx < m.Rows(); ridx++) {
                    s << m(ridx, 0);
                    for (int cidx = 1; cidx < m.Cols(); cidx++)
                        s << "\t" << m(ridx, cidx);
                    s << std::endl;
                }
                return s;
            }

        public:
            //! @brief 获取矩阵行数
            inline int Rows() const { return _numRows; }
            //! @brief 获取矩阵列数
            inline int Cols() const { return _numCols; }
            //! @brief 获取矩阵元素数量
            inline int NumDatas() const { return Rows() * Cols(); }
            //! @brief 获取矩阵总字节数
            inline int NumBytes() const { return NumDatas() * sizeof(Scalar); }
            //! @brief 获取矩阵元素字节数
            inline int DataBytes() const { return sizeof(Scalar); }
            //! @brief 获取矩阵的存储方式
            inline EStorageOptions GetStorageOption() const { return _option; }
        private:
            //! @brief 矩阵数据缓存
            std::vector<Scalar> mData;
            //! @brief 矩阵视图
            MatView mView;
    };

}


#endif
