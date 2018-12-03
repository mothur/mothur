# Mothur install instructions

Either download the precompiled binaries or compile from the source code. More detailed installation instructions are on [the mothur wiki](https://www.mothur.org/wiki/Installation).

## Download precompiled binaries

The easiest way to get mothur is to download the release from [GitHub]([GitHub](https://github.com/mothur/mothur/releases)), unzip it, and you're ready to run Mothur.

## Compile mothur from source
(For Unix-based operating systems.)

Download the mothur [source code](https://github.com/mothur/mothur).

### Compiling with Boost:

#### 1. Install dependencies.

 You will need to install the following dependencies for Boost if not already on your machine:

* bzip2
* bzip2-devel
* libz
* python
* zlib-devel

You can use a package manager such as yum, apt-get, homebrew, or conda.

#### 2. Download [Boost](http://www.boost.org).

#### 3. Follow their install [instructions]( http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html#easy-build-and-install):

```
tar -xzvf boost_versionNumber.tar.gz
cd boost_versionNumber/
./bootstrap.sh --prefix=/desired/install/path
./b2 install
```

#### 4. Compile mothur:

```
cd /path/to/mothur
make
```

If you get linking errors, it is likely because the zlib files were not found. You may need to add gzip.cpp and zlib.cpp to the source folder of mothur.  They are located in boost_versionNumber/libs/iostreams/src/.


### Compiling with HDF5:

#### 1. Download and install [HDF5]( https://portal.hdfgroup.org/display/support/HDF5+1.10.3).

```
tar -xzvf hdf5-1.10.3.tar.gz
cd hdf5-1.10.3
./configure --prefix=/desired/install/path --enable-cxx --enable-static --disable-shared
make check
make install
```

#### 2. Edit the mothur makefile.

```
cd /path/to/mothur
```

Open the makefile in your preferred text editor:

```
vi Makefile
```

And edit the HDF5 filepaths:

```
HDF5_LIBRARY_DIR ?= "/path/to/hdf5/lib"
HDF5_INCLUDE_DIR ?= "/path/to/hdf5/include"
```

Save and close the makefile. (vi command `:wq`)

#### 3. Compile mothur.
```
make
```
