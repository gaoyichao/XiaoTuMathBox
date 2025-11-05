#ifndef __XTMB_GEO_POLYGON_H__
#define __XTMB_GEO_POLYGON_H__

#include <XiaoTuDataBox/VectorList.hpp>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <vector>
#include <iostream>

namespace xiaotu {


    /**
     * @brief 给定三个点，计算三角形 a->b->c 的面积的两倍
     */
    template <typename DataType>
    DataType __Area2__(Point2<DataType> const & a,
                   Point2<DataType> const & b,
                   Point2<DataType> const & c)
    {
        return (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
    }

    /**
     * @brief 二维欧式空间下的多边形
     */
    template <typename DataType>
    class Polygon {
        public:
            using PointType = Point2<DataType>;
            using VerticesList = xiaotu::VectorList<PointType>;
            using Vertex = typename VerticesList::Node;
        
        public:
            /**
             * @brief 默认构造
             */
            Polygon() {}

            /**
             * @brief 通过一组顶点坐标构造
             */
            Polygon(std::initializer_list<PointType> && li)
            {
                Add(std::move(li));
            }

            /**
             * @brief 增加一个顶点
             */
            Vertex & Add(DataType x, DataType y, Vertex * head = nullptr)
            {
                return mVertices.Add({x, y}, head);
            }

            /**
             * @brief 添加一组顶点
             */
            Polygon & Add(std::initializer_list<PointType> && li, Vertex * head = nullptr)
            {
                for (auto rit = li.begin(); rit != li.end(); rit++)
                    mVertices.Add(*rit, head);
                return *this;
            }

            /**
             * @brief 面积的两倍
             * 
             * @param [in] idx 起始顶点
             */
            DataType Area2(int idx = 0) const
            {
                DataType sum = 0;

                Vertex const & p = mVertices[idx];
                if (!p.MoreThan2())
                    return 0;

                auto const * a = p.NextPtr();
                do {
                    sum += __Area2__(p.Data(), a->Data(), a->NextPtr()->Data());
                    a = a->NextPtr();
                } while (a->NextPtr() != &p);

                return sum;
            }

            /**
             * @brief 计算以 idx 开始的多边形面积
             */
            DataType Area(int idx = 0) const
            {
                return 0.5 * Area2(idx);
            }
            /**
             * @brief 顶点的数量
             */
            size_t NumVertices() const { return mVertices.Nodes().size(); }

            /**
             * @brief 获取指定索引的顶点
             */
            Vertex & operator [] (int idx) { return mVertices[idx]; }


        private:
            //! 多边形的顶点列表
            VerticesList mVertices;
    };

}

#endif
