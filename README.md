# ModelReader

This is a single C++ Model.h file that can be used to read in 3D models in .obj/.mtl format. It also reads in an textures that are listed in the .mtl files.

This implementation is very fast.  It can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7.

There are some limitations:
1) It does not implement all .obj features.  
2) In order to avoid dependencies on image libraries such as libjpg, it requires that all textures are first converted to .bmp format.  This can be done in one command by saying: cd texture; mogrify *.{jpg,png} -format bmp.  Model.h will auatomatically change .jpg et al. to .bmp so there is no need to edit the .mtl files.
