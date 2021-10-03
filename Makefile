
BASE_DIR = $(shell pwd)
SRC_DIR = ${BASE_DIR}/src
INC_DIR = ${BASE_DIR}/include
TEST_DIR = ${BASE_DIR}/test

BUILD_DIR = ${BASE_DIR}/build
INSTALL_CFG_FILE = ${BUILD_DIR}/all

INCLUDE_PATH += -I${INC_DIR}
INCLUDE_PATH += -I/opt/maiya/local/include/eigen3
export INCLUDE_PATH

SUBDIRS += src test 

${INSTALL_CFG_FILE}: build ${SUBDIRS}

src: FORCE
	@${MAKE} -C $@ OUTPUT_PATH=${BASE_DIR}/build

test: FORCE
	@${MAKE} -C $@ OUTPUT_PATH=${BASE_DIR}/build


build:
	mkdir build

FORCE:

install:
	-cp ${INC_DIR}/XiaoTuMathBox /usr/local/include -R
	-cp ${BUILD_DIR}/libXiaoTuMathBox.a /usr/local/lib

uninstall:
	-rm /usr/local/include/XiaoTuMathBox -r
	-rm /usr/local/lib/libXiaoTuMathBox.a

clean: clean_subdirs
	-rm build -r

-include ./build_tools/subdirs.mk


