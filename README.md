# ModelReader

This is a single C++ Model.h file that reads in 3D models in .obj/.mtl format. It also reads in any textures that are listed in the .mtl files.

The implementation uses a hand-optimized parser that can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7-6600K with one CPU core utilized.  That's fast enough to not bother saving to a binary format.  

The implementation allocates contiguous arrays (positions, normals, texcoords, vertexes, polygons, materials, textures, texels, strings) that contain no pointers, only indices into other arrays.  This means that the arrays can be copied to a GPU without editing.  They can also be stored directly to a binary file.  A hdr structure holds the lengths of the arrays.

There are also a couple std::map's (name_to_obj, name_to_tex) that take a string and look up the object/texture.

Aside from the constructor and destructor, there are no other public methods.  All arrays and maps are public.   This is a bare-bones operation here.

If the constructor encounters an error, it does not raise an exception.  Instead, it sets is_good to false and sets error_msg to a useful message.  So the caller should check is_good in the newly created Model before proceeding.

Refer to Model.h for further instructions.

There are some limitations:

1) It does not implement all .obj features.  In particular, it does not implement curves.
2) In order to avoid dependencies on image libraries such as libjpg, you must first convert all textures to .bmp format.  This can be done in one command by saying: cd textures; mogrify *.{jpg,png} -format bmp.  Model.h will automatically change .jpg et al. file extensions to .bmp so there is no need to edit the .mtl files.
3) Textures are not mipmapped - only the highest resolution is available in memory. As in the .bmp file, each row of texels is padded to a 4B boundary.  All texels have 8-bit RGB components.
4) Integers are stored as uint32_t.  Floating-point is stored as float.  However, you can change this by editing the typedefs for uint and real near the top of Model.h. 
3) .fbx is not yet supported and will be the next format to add.

Bob Alfieri
