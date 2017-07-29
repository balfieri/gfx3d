# ModelReader

This is a single C++ Model.h file that reads in 3D models in .obj/.mtl format. It also reads in an textures that are listed in the .mtl files.

This implementation is very fast.  It can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7.

Refer to the comments in Model.h for usage instructions.

There are some limitations:
1) It does not implement all .obj features.  In particular, it does not implement curves.
2) In order to avoid dependencies on image libraries such as libjpg, you must first convert all textures to .bmp format.  This can be done in one command by saying: cd textures; mogrify *.{jpg,png} -format bmp.  Model.h will automatically change .jpg et al. file extensions to .bmp so there is no need to edit the .mtl files.
