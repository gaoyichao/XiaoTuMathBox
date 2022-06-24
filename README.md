# 小土的数学工具箱 

* Ubuntu 18.04
* g++ 7.4.0
* gcc 7.4.0
* make 4.1

## 安装:

```
git clone https://github.com/gaoyichao/XiaoTuMathBox.git
cd XiaoTuMathBox
make
# 默认安装路径 /usr/local
# 如需调整，自行更改 Makefile
make install
# 卸载
make uninstall
```

## 依赖:

### 1. gtest
```
$ git clone https://github.com/google/googletest.git
$ cd gooletest
$ mkdir build
$ cd build
$ cmake ..
$ make
$ sudo make install
```

## 关联

* [几何](https://gaoyichao.com/Xiaotu/?book=几何&title=index)
* [机器人控制的数学物理基础](https://gaoyichao.com/Xiaotu/?book=math_physics_for_robotics&title=index)



