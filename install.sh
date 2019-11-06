ROOT=$PWD

CORES=1

if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform        
    CORES=$(sysctl -n hw.ncpu)
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under GNU/Linux platform
    CORES=$(nproc)
fi

#we check if the libspoa exist
if [ ! -s "${ROOT}/libs/spoa/build/lib/libspoa.a" ] 
then
	echo "Building SPOA lib"
	cd ${ROOT}/libs/spoa
	mkdir build
	cd build
	cmake ..
	make
	cd ${ROOT} 
else
	echo "SPOA lib was build before"
fi
#we compile the liger code
make -j ${CORES} 
#we check if liger was build
if [ -s "liger" ]
then 
	echo "liger binary was built"
else
	echo "problems building liger"
fi
