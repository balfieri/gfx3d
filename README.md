# ModelReader

This is a single C++ Model.h file that can be used to read in a 3D model in .obj format. 

This implementation is very fast.  It can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7.

There are some limitations:
1) It does not implement all .obj features (yet).
2) In order to avoid dependencies on image libraries, it requires that all textures are first converted to .bmp format.  This can be done in one command by saying: cd texture; mogriphy *.{jpg,png} -format bmp
