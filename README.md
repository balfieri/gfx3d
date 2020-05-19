# Model.h

This is a single C++ Model.h file (plus stb_image.h and stb_image_write.h) that reads in 3D models in .obj/.mtl format
and helps prepare them for rendering, which is handled separately.  You can think of this as the MECHANISM of 3D
graphics and rendering as the POLICY.

It also reads in any textures that are listed in the .mtl files and supports ASTC compression.
It can optionally generate mipmaps, and it can also generate a BVH tree for ray tracing.

The implementation uses a hand-optimized parser that can read in a 10M polygon model of San Miguel in less than 4 seconds on a 4GHz i7-6700K with one CPU core active.  The implementation can also write out a compressed or uncompressed file, the latter of which can be loaded
instantaneously.

The implementation allocates contiguous arrays (instances, objects, polygons, vertexes, positions, normals, texcoords, materials, textures, texels, strings) that contain no pointers, only indices into other arrays.  This means that the arrays can be copied to a GPU without editing.  They can also be stored directly to a binary file.  A hdr structure holds the lengths of the arrays.  hdr->byte_cnt is set to the total number of bytes in the header and arrays.

The Instance type allows instancing over other models, including a 4x4 Matrix transformation.

There are a couple std::map's (name_to_obj, name_to_tex) that take a string and return a pointer to the Object/Texture structure in the objects/textures array.  This are filled in only during the original parsing of the model files.

Textures are not mipmapped by default, but Model can also generate mipmaps.  The texels are stored in RGB8, RGBA8, L8, LA8, or ASTC (compressed) formats.

Errors do not raise exceptions.  Instead the constructor sets is_good to false and sets error_msg to a useful string.  So the caller should check is_good in the newly created Model before proceeding to use the object.

This has been tested on macOS and Linux.  It should work in any UNIX-like environment including Cygwin. It requires the -std=c++17 
(or later) compiler switch.

Refer to Model.h for further instructions.

This is all open-source.  Refer to the LICENSE.md for licensing details.

There are some limitations:

1) It does not implement all .obj features.  In particular, it does not implement curves or ad hoc vertex attributes.
2) Integers are stored as uint32_t except for texel offsets which are uint64_t.  Floating-point is stored as float.  
However, you can change this by editing the typedefs for uint and real near the top of Model.h.  Model.h will likely
provide ifdefs to change these rather than use templates, which can be a pain to deal with. 

Bob Alfieri<br>
Chapel Hill, NC
