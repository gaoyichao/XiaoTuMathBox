#ifndef XTMB_LA_MATRIX_COMMA_H
#define XTMB_LA_MATRIX_COMMA_H

#ifndef XTMB_LINEAR_ALGIBRA_H
#error "不要直接应用该头文件, #include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>"
#endif


namespace xiaotu::math {

    //! @brief 类似 Eigen 中的赋值语法器
    //!
    //! matrix << 1, 2,
    //!           3, 4;
    //!
    //! 参考 eigen3/Eigen/src/Core/CommaInitializer.h
    //!
    //! MatrixBase<Derived>
    template <typename Derived> 
    class MatrixComma
    {
        public:
            typedef typename Traits<Derived>::Scalar Scalar;

        public:
            MatrixComma(Derived & m, Scalar const & s)
                : target_(m), row_idx_(0), col_idx_(0), cur_row_block_size_(1)
            {
                assert(m.Cols() > 0);
                assert(m.Rows() > 0);

                target_.At(row_idx_, col_idx_) = s;
                col_idx_++;
            }

            template <typename OtherDerived>
            MatrixComma(Derived & m, MatrixBase<OtherDerived> const & s)
                : target_(m), row_idx_(0), col_idx_(0), cur_row_block_size_(s.Rows())
            {
                assert(m.Cols() > 0);
                assert(m.Rows() > 0);

                target_.Assign(row_idx_, col_idx_, s);
                col_idx_ += s.Cols();
            }

            MatrixComma & operator , (Scalar const & s)
            {
                if (target_.Cols() == col_idx_) {
                    row_idx_ += cur_row_block_size_;
                    col_idx_ = 0;
                    cur_row_block_size_ = 1;
                    assert(row_idx_ < target_.Rows());
                }
                assert(1 == cur_row_block_size_);
                assert(col_idx_ < target_.Cols());
                target_.At(row_idx_, col_idx_) = s;
                col_idx_++;

                return *this;
            }

            template <typename MatrixA>
            MatrixComma & operator , (MatrixBase<MatrixA> const & m)
            {
                if (target_.Cols() == col_idx_) {
                    row_idx_ += cur_row_block_size_;
                    col_idx_ = 0;
                    cur_row_block_size_ = m.Rows();
                    assert((row_idx_ + cur_row_block_size_) <= target_.Rows());
                }
                assert(m.Rows() == cur_row_block_size_);
                assert((col_idx_ + m.Cols()) <= target_.Cols());
                target_.Assign(row_idx_, col_idx_, m);
                col_idx_ += m.Cols();

                return *this;
            }

        private:
            //! @brief 目标矩阵
            Derived & target_;
            //! @brief 当前行更新索引
            int row_idx_;
            //! @brief 当前列更新索引
            int col_idx_;
            //! @brief 当前锁定列高度
            int cur_row_block_size_;
    };

    template <typename Derived> 
    MatrixComma<Derived> MatrixBase<Derived>::operator
    << (typename MatrixBase<Derived>::Scalar const & s)
    {
        return MatrixComma<Derived>(*static_cast<Derived*>(this), s);
    }


    template <typename Derived> 
    template <typename OtherDerived> 
    MatrixComma<Derived> MatrixBase<Derived>::operator
    << (MatrixBase<OtherDerived> const & m)
    {
        return MatrixComma<Derived>(*static_cast<Derived*>(this), m);
    }

}

#endif
