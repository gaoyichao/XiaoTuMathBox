#ifndef XTMB_LA_DECLARATIONS_H
#define XTMB_LA_DECLARATIONS_H


namespace xiaotu {
namespace math {

    template<typename T> struct traits;
    
    //! @brief 矩阵视图
    template <typename T, int numRows, int numCols, EStorageOptions option = EStorageOptions::eColMajor>
    class MatrixView;

    //! @brief 矩阵相关的基类
    template<typename Derived> class MatrixViewBase;
}
}

#endif
