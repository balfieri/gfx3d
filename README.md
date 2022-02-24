# Model.h

This is a single C++ Model.h file that reads in 3D models in .obj/.mtl format
and helps prepare them for rendering, which is handled separately.  You can think of this as the MECHANISM of 3D
graphics and rendering as the POLICY. In my proprietary path tracer, Model.h accounts for about 75% of the source lines, and I'd like to get that closer to 90% someday.

It reads in any textures that are listed in the .mtl files and supports ASTC compression.
It optionally generate mipmaps, and it can also generate a BVH tree for ray tracing.

It can read in NanoVDB volumetric databases, which is a simple format derived from OpenVDB that NVidia has open-sourced.

The implementation uses a hand-optimized parser and can also write out a single compressed or uncompressed .model file, 
the latter of which can be loaded instantaneously without translation.

The implementation allocates contiguous arrays (instances, objects, polygons, vertexes, positions, normals, texcoords, materials, textures, texels, strings, volumes, volume_grids, voxels) that contain no pointers, only indices into other arrays.  This means that the arrays can be copied to a GPU without editing and moved around.  They can also be stored directly to a binary file without translation.  A hdr structure holds the lengths of the arrays.  

Textures are not mipmapped by default, but Model can also generate mipmaps.  The texels are stored in RGB8, RGBA8, L8, LA8, or ASTC (compressed) formats.

A model may optionally instance one or more other submodels, each with a per-instance 4x4 matrix transformation.  This instancing is done within a top-level (NVidia Falcor) .fscene file which Model.h knows how to parse.  The .fscene format supports other global scene information such as sky boxes, background, ambient, tone mapping, cameras, and user-inserted light sources.  A submodel may also instance lower-level submodels.

Errors do not raise exceptions.  Instead the constructor sets is_good to false and sets error_msg to a useful string.  So the caller should check is_good in the newly created Model before proceeding to use the object.

This has been tested on macOS and Linux.  It should work in any UNIX-like environment including Cygwin. It requires the -std=c++17 
(or later) compiler switch.

Refer to Model.h for further instructions on usage.

There is an optional header called sys.h that #includes Model.h and provides some additional shorthands and utility functions which are highlighted at the top fo the file. 
If you decide to use this, you can just #include sys.h in your program rather than Model.h.

This is all open-source.  Refer to the LICENSE.md for licensing details.  We inline two well-known public-domain files, stb_image.h and
stb_image_write.h so you don't have to worry about that.  Their licensing information is listed at the top of Model.h.

There are some limitations:

1) It does not implement all .obj features.  In particular, it does not implement curves or ad hoc vertex attributes.
2) It does not support FBX files, but the graphics world seems to be standardizing on USD, so USD will get added first and we might just skip FBX.

Future features:

1) Generalize the Object class to subsume Instance and to provide a scene graph support (i.e., hierarchical objects).
2) Add a Rendering.h to do at least traditional rendering of Model.h models in Vulkan/Metal.

Bob Alfieri<br>
Chapel Hill, NC
