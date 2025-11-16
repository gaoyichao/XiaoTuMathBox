#ifndef __XTMB_GEO_POLYGON_H__
#define __XTMB_GEO_POLYGON_H__

#include <XiaoTuDataBox/VectorList.hpp>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <vector>
#include <iostream>

namespace xiaotu {


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
                    sum += xiaotu::Area2(p.Data(), a->Data(), a->NextPtr()->Data());
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
             * @brief 判定 a,b 是否为多边形的对角线
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool Diagnalie(int a, int b)
            {
                return Diagnalie(mVertices[a], mVertices[b]);
            }

            /**
             * @brief 判定 a,b 是否为多边形的对角线
             * 
             * 如果 ab 是多边形的一条边，该函数也会认为它是对角线
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool Diagnalie(Vertex const & a, Vertex const & b)
            {
                Vertex const * c = &a;
                Vertex const * c1 = nullptr;

                do {
                    c1 = c->NextPtr();
                    if ((c != &a) && (c1 != &a) && (c != &b) && (c1 != &b)
                     && (Intersect(*a, *b, **c, **c1)))
                        return false;
                    c = c1;
                } while (c != &a);
                return true;
            }
            
            /**
             * @brief 判定 a,b 是否严格地在多边形内部
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool InCone(int a, int b)
            {
                return InCone(mVertices[a], mVertices[b]);
            }

            /**
             * @brief 判定 a,b 是否严格地在多边形内部
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool InCone(Vertex const & a, Vertex const & b)
            {
                Vertex const & a0 = a.Prev();
                Vertex const & a1 = a.Next();

                if (LeftOn(*a, *a1, *a0))
                    return Left(*a, *b, *a0) && Left(*b, *a, *a1);
                return !(LeftOn(*a, *b, *a1) && LeftOn(*b, *a, *a0));
            }

            /**
             * @brief 判定 a,b 是否为多边形的对角线
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool Diagnal(int a, int b)
            {
                return Diagnal(mVertices[a], mVertices[b]);
            }

            /**
             * @brief 判定 a,b 是否为多边形的对角线
             * 
             * @param [in] a 一个顶点
             * @param [in] b 另一个顶点
             */
            bool Diagnal(Vertex const & a, Vertex const & b)
            {
                return InCone(a, b) && InCone(b, a) && Diagnalie(a, b);
            }


            /**
             * @brief 标记各个顶点是否为耳朵尖，用于计算三角分割的准备动作。
             * 
             * @param [out] earflags 对应每个顶点，true 表示是耳朵尖，false 不是
             */
            void EarInit(std::vector<bool> & earflags)
            {
                earflags.resize(NumVertices());
                earflags.assign(NumVertices(), false);

                Vertex const * head = &mVertices[0];
                Vertex const *v1 = head;
                do {
                    Vertex const *v2 = v1->NextPtr();
                    Vertex const *v0 = v1->PrevPtr();
                    earflags[v1->Index()] = Diagnal(*v0, *v2);
                    v1 = v1->NextPtr();
                } while (v1 != head);
            }


            /**
             * @brief 对多边形进行三角分割。
             * 
             * @todo 目前该函数调用之后，会改变链表的拓扑关系，后续会出一个不改变的版本
             */
            void Triangulate(std::vector<size_t> & trangulate_list)
            {
                std::vector<bool> earflags;
                EarInit(earflags);

                Vertex * head = &mVertices[0];
                int n = NumVertices();
                while (n > 3) {
                    Vertex * v2 = head;

                    do {
                        if (earflags[v2->Index()]) {
                            Vertex * v1 = v2->PrevPtr();
                            Vertex * v0 = v1->PrevPtr();
                            Vertex * v3 = v2->NextPtr();
                            Vertex * v4 = v3->NextPtr();

                            trangulate_list.push_back(v2->Index());
                            trangulate_list.push_back(v1->Index());
                            trangulate_list.push_back(v3->Index());
                            
                            earflags[v1->Index()] = Diagnal(*v0, *v3);
                            earflags[v3->Index()] = Diagnal(*v1, *v4);

                            v1->NextIdx(v3->Index());
                            v3->PrevIdx(v1->Index());
                            head = v3;
                            n--;
                            break;
                        }

                        v2 = v2->NextPtr();
                    } while (v2 != head);
                }

                assert(3 == n);
                trangulate_list.push_back(head->Index());
                trangulate_list.push_back(head->PrevIdx());
                trangulate_list.push_back(head->NextIdx());
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
