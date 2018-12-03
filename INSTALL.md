#For Unix / Linux / Centos

Compiling with Boost:
1. Boost requires some things to installed on your machine already.  Most come standard on many flavors on Unix, but you may need to install the devel packages. Install libz, bzip2 and python, if its not on your machine, including zlib-devel, bzip2-devel. 
	This can be easily done with yum or apt get.

2. Download Boost, http://www.boost.org

3. Follow their install instructions, http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html#easy-build-and-install

	./bootstrap.sh --prefix=/usr/local/
	./b2 install

4. Run make. If you get a linking errors, it is likely because the zlib files were not found correctly. You may need to add gzip.cpp and zlib.cpp to the source folder of mothur.  They are located in the boost_versionNumber/libs/iostreams/src/gzip.cpp.


Compiling with HDF5:

1. Download the tar.gz from https://portal.hdfgroup.org/display/support/HDF5+1.10.3

2. tar -xf hdf5-1.10.3.tar

3. cd hdf5-1.10.3

4. ./configure --prefix=/usr/local --enable-cxx --enable-static --disable-shared

5. make check

6. make install

7. Move LIBS libhdf5_hl_cpp.a, libhdf5_cpp.a, libhdf5_hl.a, libhdf5.a to desired location and replace {HDF5_LIBRARY_DIR} with that location in mothur's makefile. ie. /usr/local/lib HDF5_LIBRARY_DIR="/usr/local/lib"

8. Set include file location in mothur's makefile. ie. /usr/local HDF5_INCLUDE_DIR="/usr/local/"



