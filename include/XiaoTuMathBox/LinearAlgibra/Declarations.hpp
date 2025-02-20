#ifndef XTMB_LA_DECLARATIONS_H
#define XTMB_LA_DECLARATIONS_H


namespace xiaotu::math {

    //! @brief 矩阵相关的类型萃取器
    template<typename T>
    struct Traits;

    //! @brief 矩阵视图
    template <typename T, int numRows, int numCols, EAlignType align = EAlignType::eColMajor>
    class MatrixView;

    //! @brief 列向量视图
    template <typename _Scalar, int _numRows, EAlignType _align = EAlignType::eColMajor>
    using VectorView = MatrixView<_Scalar, _numRows, 1, _align>;

    //! @brief 行向量视图
    template <typename _Scalar, int _numCols, EAlignType _align = EAlignType::eColMajor>
    using RowVectorView = MatrixView<_Scalar, 1, _numCols, _align>;

    //! @brief 稠密矩阵
    template <typename T, int numRows, int numCols,
              EAlignType align = EAlignType::eColMajor,
              EStoreType store = EStoreType::eStoreVector>
    class Matrix;

    //! @brief 稠密列向量
    template <typename _Scalar, int _numRows,
                EAlignType _align = EAlignType::eColMajor,
                EStoreType _store = EStoreType::eStoreVector>
    using Vector = Matrix<_Scalar, _numRows, 1, _align, _store>;

    //! @brief 稠密行向量
    template <typename _Scalar, int _numCols,
                EAlignType _align = EAlignType::eColMajor,
                EStoreType _store = EStoreType::eStoreVector>
    using RowVector = Matrix<_Scalar, 1, _numCols, _align, _store>;

    /**
     * @brief 稠密矩阵
     * 
     * 以 std::vector<Scalar> 保存数据，适用与矩阵尺寸较大的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用小。
     * 矩阵尺寸小时，大量使用，将导致内存碎片化。
     */
    template <typename T, int numRows, int numCols,
            EAlignType align = EAlignType::eColMajor>
    using VMatrix = Matrix<T, numRows, numCols, align, EStoreType::eStoreVector>;

    /**
     * @brief 稠密矩阵
     * 
     * 以 Scalar[] 保存数据，适用与矩阵尺寸较小的情况。
     * 矩阵尺寸大时，用作局部变量，栈空间占用大。
     */
    template <typename T, int numRows, int numCols,
            EAlignType align = EAlignType::eColMajor>
    using AMatrix = Matrix<T, numRows, numCols, align, EStoreType::eStoreArray>;
    
}

#endif
