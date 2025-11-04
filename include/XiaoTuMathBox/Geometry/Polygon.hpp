#ifndef __XTMB_GEO_POLYGON_H__
#define __XTMB_GEO_POLYGON_H__

#include <XiaoTuDataBox/VectorList.hpp>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <vector>
#include <iostream>

namespace xiaotu::math {

    /**
     * @brief 二维欧式空间下的多边形
     */
    template <typename DataType>
    class Polygon {
        public:
            using VerticesList = xiaotu::data::VectorList<Point2<DataType>>;
            using Vertex = typename VerticesList::Node;
        
        public:
            /**
             * @brief 默认构造
             */
            Polygon() {}

            /**
             * @brief 顶点的数量
             */
            size_t NumVertices() const { return mVertices.Nodes().size(); }

            /**
             * @brief 增加一个顶点
             */
            Vertex & Add(DataType x, DataType y, Vertex * head = nullptr)
            {
                return mVertices.Add({x, y}, head);
            }

            Vertex & operator [] (int idx) { return mVertices[idx]; }

        private:
            //! 多边形的顶点列表
            VerticesList mVertices;
    };

}

#endif
