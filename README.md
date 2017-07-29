# ModelReader

This is a single C++ Model.h file that reads in 3D models in .obj/.mtl format. It also reads in any textures that are listed in the .mtl files.

The implementation uses a hand-optimized parser that can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7-6600K with one CPU core utilized.  That's fast enough to not bother with saving a to binary format.  However, I might add multi-threading at some point.

The implementation allocates contiguous arrays (positions, normals, texcoords, vertexes, polygons, materials, textures, strings) that contain no pointers, only index into other arrays.  This means that the arrays can be copied into a GPU or written out to a binary file without having to edit them.

Refer to the comments at the top of Model.h for usage instructions.

There are some limitations:
1) It does not implement all .obj features.  In particular, it does not implement curves.
2) In order to avoid dependencies on image libraries such as libjpg, you must first convert all textures to .bmp format.  This can be done in one command by saying: cd textures; mogrify *.{jpg,png} -format bmp.  Model.h will automatically change .jpg et al. file extensions to .bmp so there is no need to edit the .mtl files.
3) .fbx is not yet supported and will be the next format to add.

Bob Alfieri
