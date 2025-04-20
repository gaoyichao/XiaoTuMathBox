# 小土的数学工具箱 

* Ubuntu 18.04
* g++ 7.4.0
* gcc 7.4.0
* make 4.1

## 安装:

```
git clone https://github.com/gaoyichao/XiaoTuMathBox.git
cd XiaoTuMathBox
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/local ..
make
# 安装
make install
# 卸载
make uninstall
# 单元测试
make test
```

## 依赖:

### 1. gtest
```
git clone https://github.com/google/googletest.git
cd gooletest
mkdir build
cd build
cmake ..
make
sudo make install
```

### 2. eigen
```
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cd build
cmake ..
make
sudo make install
```

### 3. glog
```
git clone https://github.com/google/glog.git
cd glog
mkdir build
cd build
cmake ..
make
sudo make install
```

### 4. gflags
```
git clone https://github.com/gflags/gflags.git
cd gflags/
mkdir build
cd build/
cmake ..
make
sudo make install
```


## 关联

* [线性代数](http://gaoyichao.com/Xiaotu/?book=algebra&title=index)
* [几何](https://gaoyichao.com/Xiaotu/?book=几何&title=index)
* [机器人控制的数学物理基础](https://gaoyichao.com/Xiaotu/?book=math_physics_for_robotics&title=index)



