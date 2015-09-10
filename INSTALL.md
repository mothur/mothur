#For Unix / Linux / Centos

1. Boost requires some things to installed on your machine already.  Most come standard on many flavors on Unix, but you may need to install the devel packages. Install libz, bzip2 and python, if its not on your machine, including zlib-devel, bzip2-devel. 
	This can be easily done with yum or apt get.

2. Download Boost, http://www.boost.org

3. Follow their install instructions, http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html#easy-build-and-install

	./bootstrap.sh --prefix=/usr/local/
	./b2 install

4. Run make. If you get a linking errors, it is likely because the zlib files were not found correctly. You may need to add gzip.cpp and zlib.cpp to the source folder of mothur.  They are located in the boost_versionNumber/libs/iostreams/src/gzip.cpp.


