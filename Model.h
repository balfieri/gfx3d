// Copyright (c) 2017-2021 Robert A. Alfieri
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
//
// Model.h - 3D model loading with one .h file
//
// We inline these two well-known public-domain headers at the end of this file 
// so you still need only include Model.h (and you should not include those separately beforehand):
//
// https://github.com/nothings/stb (MIT license)
//
// stb_image.h     - v2.06 - public domain image loader 
//                           no warranty implied; use at your own risk
//
// stb_image_write - v0.98 - public domain image writer
//                           writes out PNG/BMP/TGA images to C stdio - Sean Barrett 2010
//                           no warranty implied; use at your own risk
//
// Copyright (c) 2017 Sean Barrett
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// How to Model.h:
//
//     1) #include "Model.h"
//
//        Model * model = new Model( "~/models/sanmiguel/sanmiguel.obj" );   
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//        model->write( "~models/sanmiguel/sanmiguel.model" );   // will write out the self-contained binary model
//                                                               // default is compressed, add "false" as 2nd arg for uncompressed
//
//     2) After that, you can quickly read in the single binary model file using:
//
//        Model * model = new Model( "~models/sanmiguel/sanmiguel.model", false );  // false means uncompressed
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//
//     3) If you want Model to generate mipmap textures, add Model::MIPMAP_FILTER::BOX as a second argument 
//        to the Model() constructor in (1) to get a box filter.  Other filters may be added later.
//        The texture for mip level 0 is the original texture.  The texels for the other mip levels
//        follow immediately with no padding in between.  The original width and height need not
//        be powers-of-2 or equal to each other.  Each mip level has 
//        width=min(1, prev_width>>1) x height=min(1, prev_height>>1) texels.
//        The last mip level will always contain 1x1 texels.  Model does not currently store
//        the number of texels in each level, so you'll need to compute those on the fly as you
//        try to obtain the starting offset for a given level.
//
//     4) If you want Model to generate compressed textures for you, it will use ARM's astcenc program to encode
//        the textures using ASTC which is a modern, high-quality compression format that is often 10-20x smaller than uncompressed.
//        You should pass in the original (ideally uncompressed) file name, e.g., image1.png, so that we can encode using the best
//        version of the image.  Model will create files in the same directory as image1.png named 
//        image1.astc for mip level 0 (finest), image1.1.astc for mip level 1, etc.  This way, if these
//        files already exist, Model will not need to re-generate them, which can otherwise take quite a bit of time.
//        Your .mtl files may also refer to .astc files for mip 0, which can be named image1.astc etc.
//        In this case, Model will generate image1.1.astc, etc. for the other mip levels if they do not
//        already exist.  It is better to let Model generate all mip levels from the original image file OR
//        generate them all yourself so they are ready to go.  In practice, I have not confirmed that
//        there's any noticeable difference.
//
//        Currently, Model uses a heuristic to pick the block size for each mip level.
//        If you would like different block sizes for different textures or mip levels, then we can add
//        an option or, for now, you can do your own texture filtering outside of Model and just be sure to name 
//        your files image1.astc, image1.1.astc, etc. so that Model will see them and not regenerate them.
//
//        As mentioned in the README.md, ASTC encode relies on the "astcenc" program being installed on your
//        system and visible on your PATH. This file does its own ASTC decode and does not rely on any external
//        programs (e.g., astcdec) or any external source code.
//
//     5) If you want Model to generate a BVH tree for you, add Model::BVH_TREE::BINARY as the third argument
//        to the Model() constructor in (1) to get a binary BVH tree.  QUAD and OCT trees are not
//        currently supported.
//
//        Note: BVH building reorders the polygons, so any polygon offsets in an Object structure will be garbage.
//              Most apps don't use the objects[] so this is not normally a problem. In any case, this will be fixed soon
//              in prepartion for hierarchical scene graphs which will make a lot of use of objects.
//
// How it works:
//
//     1) Allocate large virtual memory 1D arrays for materials, texels, positions, normals, vertexes, polygons.
//        These are allocated on a page boundary to make uncompressed writes faster.
//        These arrays are dynamically resized.
//     2) Read entire .obj file into memory (uses mmap()).
//     3) Parse .obj file using custom parser that goes character-by-character and does its own number conversions. 
//     4) Add elements to 1D arrays.  Load any .mtl file encounted in .obj file.  
//     5) Optionally generate mipmap textures.  
//     6) Optionally compress textures using ARM's ASTC compressor program: astcenc.
//     7) Optionally generate BVH tree.
//     8) Write to  uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//     9) Read from uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//        (and uses mmap() system call).
//
#ifndef _Model_h
#define _Model_h

#include <cstdint>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <mutex>

#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>
#include <float.h>

// forward decls for stb_image (some used by Model methods)
typedef unsigned char stbi_uc;
#define STBIDEF 
STBIDEF stbi_uc *stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
STBIDEF int      stbi_write_png(char const *filename, int w, int h, int comp, const void  *data, int stride_in_bytes);
STBIDEF int      stbi_write_bmp(char const *filename, int w, int h, int comp, const void  *data);
STBIDEF int      stbi_write_tga(char const *filename, int w, int h, int comp, const void  *data);
STBIDEF int      stbi_write_hdr(char const *filename, int w, int h, int comp, const float *data);
STBIDEF void     stbi_image_free(void *retval_from_stbi_load);  
STBIDEF const char *stbi_failure_reason(void);

// default scalar types used throughout
#ifndef MODEL_INT_TYPE
#define MODEL_INT_TYPE int32_t
#endif
#ifndef MODEL_INT64_TYPE
#define MODEL_INT64_TYPE int64_t
#endif
#ifndef MODEL_UINT_TYPE
#define MODEL_UINT_TYPE uint32_t
#endif
#ifndef MODEL_UINT64_TYPE
#define MODEL_UINT64_TYPE uint64_t
#endif
#ifndef MODEL_REAL_TYPE
#define MODEL_REAL_TYPE float
#define MODEL_REAL_MAX FLT_MAX
#define MODEL_REAL_EXP_W 7
#define MODEL_REAL_FRAC_W 23
#endif
#ifndef MODEL_REAL32_TYPE
#define MODEL_REAL32_TYPE float
#define MODEL_REAL32_MAX FLT_MAX
#define MODEL_REAL32_EXP_W 7
#define MODEL_REAL32_FRAC_W 23
#endif
#ifndef MODEL_REAL64_TYPE
#define MODEL_REAL64_TYPE double
#define MODEL_REAL64_MAX DBL_MAX
#define MODEL_REAL64_EXP_W 11
#define MODEL_REAL64_FRAC_W 53
#endif

// user-definable RNG function used in a few places
// this should return [0,1)
// default is good-old rand().
#ifndef MODEL_UNIFORM_FN
#define MODEL_UNIFORM_FN _model_uniform
static inline double _model_uniform() { return drand48(); }
#endif

class Model
{
public:
    typedef MODEL_INT_TYPE    _int;                  
    typedef MODEL_INT64_TYPE  _int64;                  
    typedef MODEL_UINT_TYPE   uint;                 
    typedef MODEL_UINT64_TYPE uint64;              
    typedef MODEL_REAL_TYPE   real;
    typedef MODEL_REAL32_TYPE real32;
    typedef MODEL_REAL64_TYPE real64;

    enum class MIPMAP_FILTER
    {
        NONE,                               // do not generate mipmap levels
        BOX,                                // generate mipmap levels using simple box filter (average)
    };

    enum class TEXTURE_COMPRESSION          // texture compression
    {
        NONE,                               // uncompressed
        ASTC,                               // ASTC compression 
    };

    enum class BVH_TREE
    {
        NONE,                               // do not generate BVH tree
        BINARY,                             // generate binary BVH tree
    };

    Model( std::string          top_file="",
           MIPMAP_FILTER        mipmap_filter=MIPMAP_FILTER::NONE, 
           TEXTURE_COMPRESSION  texture_compression=TEXTURE_COMPRESSION::NONE,
           BVH_TREE             bvh_tree=BVH_TREE::NONE, 
           bool                 resolve_models=true,
           bool                 deduplicate=false );
    Model( std::string model_file, bool is_compressed );
    ~Model(); 

    bool write( std::string file_path, bool is_compressed ); 
    bool replace_materials( std::string mtl_file_path );

    // start of main structures
    static const uint VERSION = 0xB0BA1f0a; // current version 

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    class real4
    {
    public:
        real c[4];
        
        inline real4( void )                               { c[0] = 0;  c[1] = 0;  c[2] = 0;  c[3] = 0;  }
        inline real4( real c0, real c1, real c2, real c3 ) { c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3; }

        inline real x() const { return c[0]; }
        inline real y() const { return c[1]; }
        inline real z() const { return c[2]; }
        inline real w() const { return c[3]; }
        inline real r() const { return c[0]; }
        inline real g() const { return c[1]; }
        inline real b() const { return c[2]; }
        inline real a() const { return c[3]; }
        
        inline const real4& operator+() const           { return *this; }
        inline real4    operator-() const               { return real4(-c[0], -c[1], -c[2], -c[3]); }
        inline real     operator[](int i) const         { return c[i]; }
        inline real&    operator[](int i)               { return c[i]; };
        inline bool     operator == ( const real4 &v2 ) const;
        inline real4&   operator += ( const real4 &v2 );
        inline real4&   operator -= ( const real4 &v2 );
        inline real4&   operator *= ( const real4 &v2 );
        inline real4&   operator *= ( const real s );
        inline real4&   operator /= ( const real4 &v2 );
        inline real4&   operator /= ( const real s );
    
        inline real     dot( const real4 &v2 ) const;
        inline real     length( void ) const;
        inline real     length_sqr( void ) const ;
        inline real     distance( const real4& other ) const;
        inline real     distance_sqr( const real4& other ) const;
        inline real4&   normalize( void );
        inline real4    normalized( void ) const;
        inline void     clamp( real min=0.0, real max=1.0 );
        inline real4    clamped( void ) const;
        inline real4    mulby2( void ) const;
        inline real4    mulby4( void ) const;
        inline real4    divby2( void ) const;
        inline real4    divby4( void ) const;
    };

    class real3d        // double-precision
    {
    public:
        real64 c[3];

        inline real3d( void )                            { c[0] = 0;  c[1] = 0;  c[2] = 0;  }
        inline real3d( real64 c0, real64 c1, real64 c2 ) { c[0] = c0; c[1] = c1; c[2] = c2; }

        inline real64 x() const { return c[0]; }
        inline real64 y() const { return c[1]; }
        inline real64 z() const { return c[2]; }
        inline real64 r() const { return c[0]; }
        inline real64 g() const { return c[1]; }
        inline real64 b() const { return c[2]; }
        
        inline const real3d& operator+() const   { return *this; }
        inline real3d   operator-() const       { return real3d(-c[0], -c[1], -c[2]); }
        inline real64   operator[](int i) const { return c[i]; }
        inline real64&  operator[](int i)       { return c[i]; };
        inline bool     operator == ( const real3d &v2 ) const;
        inline real3d&  operator += ( const real3d &v2 );
        inline real3d&  operator -= ( const real3d &v2 );
        inline real3d&  operator *= ( const real3d &v2 );
        inline real3d&  operator *= ( const real64 s );
        inline real3d&  operator /= ( const real3d &v2 );
        inline real3d&  operator /= ( const real64 s );
    
        inline real64   dot( const real3d &v2 ) const;
        inline real3d   cross( const real3d &v2 ) const;
        inline real64   length( void ) const;
        inline real64   length_sqr( void ) const ;
        inline real64   distance( const real3d& other ) const;
        inline real64   distance_sqr( const real3d& other ) const;
        inline real3d&  normalize( void );
        inline real3d   normalized( void ) const;
        inline void     clamp( real64 min=0.0, real64 max=1.0 );
        inline real3d   clamped( void ) const;
        inline real3d   mulby2( void ) const;
        inline real3d   mulby4( void ) const;
        inline real3d   divby2( void ) const;
        inline real3d   divby4( void ) const;
        inline real64   luminance( void ) const { return c[0]*0.299 +  c[1]*0.587  + c[2]*0.114; }
        inline real64   average( void ) const   { return 0.333333333333333333333333333*(c[0] + c[1] + c[2]); }
        inline real3d   sqrt( void ) const      { return real3d( std::sqrt(c[0]), std::sqrt(c[1]), std::sqrt(c[2]) ); }
        inline real3d   pow( real n ) const     { return real3d( std::pow(c[0], n), std::pow(c[1], n), std::pow(c[2], n) ); }
        inline bool     is_zero( void ) const   { return c[0] == 0.0 && c[1] == 0.0 && c[2] == 0.0; }
        inline bool     is_one( void ) const    { return c[0] == 1.0 && c[1] == 1.0 && c[2] == 1.0; }
    };

    class real3         // default precision 
    {
    public:
        real c[3];
        
        inline real3( void )                      { c[0] = 0;  c[1] = 0;  c[2] = 0;  }
        inline real3( real c0, real c1, real c2 ) { c[0] = c0; c[1] = c1; c[2] = c2; }
        inline real3( const real3d& v )           { c[0] = v.c[0]; c[1] = v.c[1]; c[2] = v.c[2]; }
        inline real3( const real4& v4 )           { c[0] = v4.c[0]; c[1] = v4.c[1]; c[2] = v4.c[2]; }

        inline real x() const { return c[0]; }
        inline real y() const { return c[1]; }
        inline real z() const { return c[2]; }
        inline real r() const { return c[0]; }
        inline real g() const { return c[1]; }
        inline real b() const { return c[2]; }
        
        inline const real3& operator+() const   { return *this; }
        inline real3    operator-() const       { return real3(-c[0], -c[1], -c[2]); }
        inline real     operator[](int i) const { return c[i]; }
        inline real&    operator[](int i)       { return c[i]; };
        inline bool     operator == ( const real3 &v2 ) const;
        inline real3&   operator += ( const real3 &v2 );
        inline real3&   operator -= ( const real3 &v2 );
        inline real3&   operator *= ( const real3 &v2 );
        inline real3&   operator *= ( const real s );
        inline real3&   operator /= ( const real3 &v2 );
        inline real3&   operator /= ( const real s );
        inline real3d   to_real3d( void ) const { return real3d(c[0], c[1], c[2]); }
    
        inline real     dot( const real3 &v2 ) const;
        inline real3    cross( const real3 &v2 ) const;
        inline real     length( void ) const;
        inline real     length_sqr( void ) const;
        inline real     distance( const real3& other ) const;
        inline real     distance_sqr( const real3& other ) const;
        inline real3&   normalize( void );
        inline real3    normalized( void ) const;
        inline void     clamp( real min=0.0, real max=1.0 );
        inline real3    clamped( void ) const;
        inline real3    mulby2( void ) const;
        inline real3    mulby4( void ) const;
        inline real3    divby2( void ) const;
        inline real3    divby4( void ) const;
        inline real     luminance( void ) const { return c[0]*0.299 +  c[1]*0.587  + c[2]*0.114; }
        inline real     average( void ) const   { return 0.333333333333333333333333333*(c[0] + c[1] + c[2]); }
        inline real3    sqrt( void ) const      { return real3( std::sqrt(c[0]), std::sqrt(c[1]), std::sqrt(c[2]) ); }
        inline real3    pow( real n ) const     { return real3( std::pow(c[0], n), std::pow(c[1], n), std::pow(c[2], n) ); }
        inline bool     is_zero( void ) const   { return c[0] == 0.0 && c[1] == 0.0 && c[2] == 0.0; }
        inline bool     is_one( void ) const    { return c[0] == 1.0 && c[1] == 1.0 && c[2] == 1.0; }
    };

    class real2
    {
    public:
        real c[2];

        inline real2( void )             { c[0] = 0;  c[1] = 0;  }
        inline real2( real c0, real c1 ) { c[0] = c0; c[1] = c1; }

        inline real x() const { return c[0]; }
        inline real y() const { return c[1]; }
        inline real u() const { return c[0]; }
        inline real v() const { return c[1]; }
        
        inline const real2& operator+() const   { return *this; }
        inline real2    operator-() const       { return real2(-c[0], -c[1]); }
        inline real     operator[](int i) const { return c[i]; }
        inline real&    operator[](int i)       { return c[i]; };
        inline real2&   operator += ( const real2 &v2 );
        inline real2&   operator -= ( const real2 &v2 );
        inline real2&   operator *= ( const real2 &v2 );
        inline real2&   operator *= ( const real s );
        inline real2&   operator /= ( const real2 &v2 );
        inline real2&   operator /= ( const real s );

        inline real     dot( const real2 &v2 ) const;
        inline real     length( void ) const;
        inline real     length_sqr( void ) const ;
        inline real     distance( const real2& other ) const;
        inline real     distance_sqr( const real2& other ) const;
        inline real2&   normalize( void );
        inline real2    normalized( void ) const;
        inline void     clamp( real min=0.0, real max=1.0 );
        inline real2    clamped( void ) const;
    };

    class Header                            // header (of future binary file)
    {
    public:
        uint        version;                // version
        uint64      byte_cnt;               // total in-memory bytes including this header
        uint64      char_cnt;               // in strings array

        uint64      obj_cnt;                // in objects array  
        uint64      poly_cnt;               // in polygons array 
        uint64      emissive_poly_cnt;      // in emissive_polygons[] array
        uint64      vtx_cnt;                // in vertexes array
        uint64      pos_cnt;                // in positions array
        uint64      norm_cnt;               // in normals array
        uint64      texcoord_cnt;           // in texcoords array

        MIPMAP_FILTER mipmap_filter;        // if NONE, there are no mip levels beyond level 0
        uint64      mtl_cnt;                // in materials array
        uint64      tex_cnt;                // in textures array
        uint64      texel_cnt;              // in texels array  

        uint64      graph_node_cnt;         // in graph_nodes array
        uint64      graph_root_i;           // index of root Graph_Node in graph_nodes array

        uint64      matrix_cnt;             // in matrixes array
        uint64      inst_cnt;               // in instances array

        uint64      bvh_node_cnt;           // in bvh_nodes array
        uint64      bvh_root_i;             // index of root BVH_Node in bvh_nodes array

        uint64      light_cnt;              // in lights array
        real        lighting_scale;         // not sure what this is yet (default: 1.0)
        real3       ambient_intensity;      // ambient light intensity (default: [0.1, 0.1, 0.1]
        real3       background_color;       // default is black, though irrelevant if there's a sky box
        bool        tex_specular_is_orm;    // specular texture components mean: occlusion, roughness, metalness
        uint        sky_box_tex_i;          // index in textures array of sky box texture
        real        env_map_intensity_scale;// name says it all
        TEXTURE_COMPRESSION texture_compression; // NONE or ASTC used (individual textures are also marked)
        real        opacity_scale;          // not sure what this is for
        real        shadow_caster_count;    // not sure what this is for
        uint        direct_V;               // number of virtual point lights (VPLs) for direct lighting reservoir sampling (default: 0 == choose)
        uint        direct_M;               // number of M candidate samples for direct lighting reservoir sampling (default: 0 == choose)
        uint        direct_N;               // number of N used samples for direct lighting reservoir sampling (default: 0 == choose)
        real        tone_white;             // tone mapping white parameter
        real        tone_key;               // tone mapping key parameter
        real        tone_avg_luminance;     // tone mapping average luminance override (if supported by application)

        uint64      camera_cnt;             // in cameras array
        uint64      initial_camera_i;       // initial active camera index (default: 0)
        uint64      frame_cnt;              // in frames array 
        uint64      animation_cnt;          // in animations array
        real        animation_speed;        // divide frame time by this to get real time (default: 1.0)

        uint64      volume_cnt;             // in volumes[] array
        uint64      volume_grid_cnt;        // in volume_grids[] array
        uint64      voxel_cnt;              // in voxels[] array  

        bool        force_tone_none;        // force using no tone-mapping?
        bool        force_tone_avg_luminance; // force the use of tone_avg_luminance above?
        bool        unused1[6];             // leave room for temporary hacks
        real        sky_radius;             // radius of background or sky box
        real        unused2[1];             // leave room for temporary hacks
        uint64_t    unused3[3];             // leave room for temporary hacks
    };

    class Object
    {
    public:
        uint        name_i;                 // index of object name in strings array
        uint        poly_cnt;               // number of polygons in this object
        uint        poly_i;                 // index of first polygon in polygons array
    };

    class HitRecord
    {
    public:
        // input state that may affect hit decisions:
        real3           F0;                 // incoming F0 (valid if F0[0] >= 0) used to decide when to stop in a volume
                                            // if you have only IOR (index of refraction), construct F0 as follows:
                                            //     a = (IOR-1)^2 / (IOR+1)^2)
                                            //     F0 = (a, a, a)

        // fields common to all types of hits
        const Model *   model;              // model of owning BVH
        real            t;                  // t parameter to get to p from origin
        real3           p;                  // surface hit point

        // it's a Polygon hit when poly_i != -1
        uint            poly_i;                 
        real            u;
        real            v;
        real            frac_uv_cov;        // UV footprint of hit surface (used by mip_level calculation)
        real3           normal;
        real3           shading_normal;
        real3           tangent;
        real3           tangent_normalized;
        real3           bitangent;

        // It's a VolumeGrid hit when volume_i != -1.
        // Currently, all grids must have the same bounding boxes.  This is so that
        // we can return a Volume when a hit() occurs for whatever reason.
        // We also return the grid_i that caused the hit() to stop.
        // Typically, the caller will want to interrogate all grids during a hit.
        //
        uint            volume_i;               
        uint            grid_i;             // grid index that caused hit to stop
        _int            voxel_xyz[3];       // voxel coordinates 
        struct
        {
            real        r;
            real64      r64;
            real3       r3;
            real3d      r3d;
            _int        i;
            _int64      i64;
        }               voxel_value;        // based on voxel_type of grid_i
    };

    class AABB                              // axis aligned bounding box
    {
    public:
        real3           _min;               // bounding box min
        real3           _max;               // bounding box max

        inline AABB( void ) {}
        inline AABB( const real3& p );             // init with one point
        inline AABB( const real3& p0, const real3& p1 );
        inline AABB( const real3& p0, const real3& p1, const real3& p2 );
        AABB( const AABB& aabb1, const AABB& aabb2 ) { *this = aabb1; expand( aabb2 ); }

        inline real3 min() const { return _min; }
        inline real3 max() const { return _max; }

        void pad( real p );
        void expand( const AABB& other );
        void expand( const real3& p );
        bool encloses( real x, real y, real z ) const { real3 p( x, y, z ); return encloses( p ); }
        inline bool is_empty( void ) const      { return _max[0] <= _min[0] || _max[1] <= _min[1] || _max[2] <= _min[2]; }
        inline real volume( void ) const        { return (_max[2]-_min[2])*(_max[1]-_min[1])*(_max[0]-_min[0]); }
        inline real3 center( void ) const       { return real3( (_min[0]+_max[0]) * 0.5, (_min[1]+_max[1]) * 0.5, (_min[2]+_max[2]) * 0.5 ); }
        inline real3 size( void ) const         { return real3( _max[0]-_min[0], _max[1]-_min[1], _max[2]-_min[2] ); }
        inline real size_max( void ) const      { real3 s = size(); return std::max( std::max( s[0], s[1] ), s[2] ); }
        inline real size_max_rcp( void ) const  { return 1.0f / size_max(); }
        inline real3 half_diag( void ) const    { return size().divby2(); }
        inline real half_area( void ) const     { real3 s = size(); return s[0]*s[1] + s[1]*s[2] + s[2]*s[0]; }
        inline real area( void ) const          { return half_area() * 2.0f; }
        bool encloses( const real3& p ) const;
        bool encloses( const AABB& other ) const;
        bool overlaps( const AABB& other ) const;
        real3 transform_relative( const real3& v ) const;
        AABB transform_relative( const AABB& other ) const;
        bool hit( const real3& origin, const real3& direction, const real3& direction_inv, real& tmin, real& tmax ) const; 
    };

    class AABBD                              // axis aligned bounding box
    {
    public:
        real3d           _min;               // bounding box min
        real3d           _max;               // bounding box max

        inline AABBD( void ) {}
        AABBD( const real3d& p );            // init with one point
        AABBD( const real3d& p0, const real3d& p1 );
        AABBD( const real3d& p0, const real3d& p1, const real3d& p2 );
        AABBD( const AABB& aabb );
        AABBD( const AABBD& aabb1, const AABBD& aabb2 ) { *this = aabb1; expand( aabb2 ); }

        inline AABB to_aabb( void ) const        { AABB aabb; aabb._min[0] = _min[0]; aabb._min[1] = _min[1]; aabb._min[2] = _min[2]; aabb._max[0] = _max[0]; aabb._max[1] = _max[1]; aabb._max[1] = _max[1]; return aabb; }

        inline real3d min() const { return _min; }
        inline real3d max() const { return _max; }

        void pad( real64 p );
        void expand( const AABBD& other );
        void expand( const real3d& p );
        inline bool   is_empty(void) const       { return _max[0] <= _min[0] || _max[1] <= _min[1] || _max[2] <= _min[2]; }
        inline real64 volume(void) const         { return (_max[2]-_min[2])*(_max[1]-_min[1])*(_max[0]-_min[0]); }
        inline real3d center( void ) const       { return real3d( (_min[0]+_max[0]) * 0.5, (_min[1]+_max[1]) * 0.5, (_min[2]+_max[2]) * 0.5 ); }
        inline real3d size( void ) const         { return real3d( _max[0]-_min[0], _max[1]-_min[1], _max[2]-_min[2] ); }
        inline real64 size_max( void ) const     { real3d s = size(); return std::max( std::max( s[0], s[1] ), s[2] ); }
        inline real64 size_max_rcp( void ) const { return 1.0 / size_max(); }
        inline real3d half_diag( void ) const    { return size().divby2(); }
        inline real64 half_area( void ) const    { real3d s = size(); return s[0]*s[1] + s[1]*s[2] + s[2]*s[0]; }
        inline real64 area( void ) const         { return half_area() * 2.0; }
        bool encloses( const real3d& p ) const;
        bool encloses( real64 x, real64 y, real64 z ) const { real3d p( x, y, z ); return encloses( p ); }
        bool encloses( const AABBD& other ) const;
        bool overlaps( const AABBD& other ) const;
        bool overlaps_triangle( const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &v2, Model::real3d * p_ptr=nullptr ) const;
        real3d transform_relative( const real3d& v ) const;
        AABBD transform_relative( const AABBD& other ) const;
        bool hit( const real3d& origin, const real3d& direction, const real3d& direction_inv, real64& tmin, real64& tmax ) const; 
    };

    class AABBI                             // axis aligned bounding box with integers
    {
    public:
        _int            _min[3];            // bounding box min
        _int            _max[3];            // bounding box max

        inline AABBI( void ) {}
        inline AABBI( const _int p[] );             // init with one point
        inline AABBI( const _int p0[], const _int p1[] );
        inline AABBI( const _int p0[], const _int p1[], const _int p2[] );
        void expand( const _int p[] );
        inline bool   is_empty(void) const       { return (_max[2]-_min[2]) <= 0 || (_max[1]-_min[1]) <= 0 || (_max[0]-_min[0]) <= 0; }
        inline uint64 volume(void) const         { return uint64(_max[2]-_min[2]+1)*uint64(_max[1]-_min[1]+1)*uint64(_max[0]-_min[0]+1); }
        bool encloses( const _int p[] ) const;
        bool encloses( _int x, _int y, _int z ) const;
    };

    class AABBU64                           // axis aligned bounding box with uint64
    {
    public:
        uint64          _min[3];            // bounding box min
        uint64          _max[3];            // bounding box max

        inline AABBU64( void ) {}
        AABBU64( const uint64 p[] );             // init with one point
        AABBU64( const uint64 p0[], const uint64 p1[] );
        AABBU64( const uint64 p0[], const uint64 p1[], const uint64 p2[] );
        void expand( const uint64 p[] );
        inline bool   is_empty(void) const       { return _max[0] <= _min[0] || _max[1] <= _min[1] || _max[2] <= _min[2]; }
        inline uint64 volume(void) const         { return (_max[2]-_min[2]+1)*(_max[1]-_min[1]+1)*(_max[0]-_min[0]+1); }
        bool encloses( const uint64 p[] ) const;
        bool encloses( uint64 x, uint64 y, uint64 z ) const;
    };

    enum class RAY_KIND 
    {
        DIRECT,
        REFLECTED,
        REFRACTED
    };

    class Ray
    {
    public: 
        inline Ray() { top = -1;}
        inline Ray( const real3& a, const real3& b, RAY_KIND kind, real32 solid_angle=0, real32 cone_angle=0 );
        inline void init_normalized( const real3& a, const real3& b, RAY_KIND kind, real32 solid_angle=0, real32 cone_angle=0 ); 

        inline const real3& origin() const            { return A; }
        inline const real3& direction() const         { return B; }
        inline const real3& direction_inv() const     { return B_inv; }
        inline const real3& normalized_direction() const { return B_norm; }
        inline RAY_KIND kind() const                  { return KIND; }
        inline real32 solid_angle() const             { return SOLID_ANGLE; }
        inline real3 point_at_parameter(real t) const;
        inline real3 normalized_point_at_parameter(real t) const;
        inline real3 endpoint() const                 { return point_at_parameter(1.0); }
        inline real32 length() const                  { return direction().length(); }

        inline real32 cone_angle() const              { return CONE_ANGLE; }
        inline real32 cone_radius(real t=1.0) const   { return (CONE_ANGLE <= 0.0 || t <= 0.0) ? 0.0 : (t*length() * std::tan(CONE_ANGLE)); }
        inline real32 cone_radius_sqr(real t=1.0) const { real32 r = cone_radius(t); return r*r; }
        inline const real3 cone_base_dir_x( void ) const;
        inline const real3 cone_base_dir_y( void ) const { return -direction(); }  // from endpoint() to origin()

        real3 A;
        real3 B;
        real3 B_inv;
        real3 B_norm;
        RAY_KIND KIND;
        real32 SOLID_ANGLE;
        real32 CONE_ANGLE;  
        int top;
    };

    class Ray64
    {
    public: 
        inline Ray64() { top = -1;}
        inline Ray64( const real3d& a, const real3d& b, RAY_KIND kind, real64 solid_angle=0, real64 cone_angle=0 );
        inline void init_normalized( const real3d& a, const real3d& b, RAY_KIND kind, real64 solid_angle=0, real64 cone_angle=0 ); 

        inline const real3d& origin() const            { return A; }
        inline const real3d& direction() const         { return B; }
        inline const real3d& direction_inv() const     { return B_inv; }
        inline const real3d& normalized_direction() const { return B_norm; }
        inline RAY_KIND kind() const                   { return KIND; }
        inline real64 solid_angle() const              { return SOLID_ANGLE; }
        inline real3d point_at_parameter(real64 t) const;
        inline real3d normalized_point_at_parameter(real64 t) const;
        inline real3d endpoint() const                 { return point_at_parameter(1.0); }
        inline real64 length() const                   { return direction().length(); }

        inline real64 cone_angle() const               { return CONE_ANGLE; }
        inline real64 cone_radius(real t=1.0) const    { return (CONE_ANGLE <= 0.0 || t <= 0.0) ? 0.0 : (t*length() * std::tan(CONE_ANGLE)); }
        inline real64 cone_radius_sqr(real t=1.0) const { real64 r = cone_radius(t); return r*r; }
        inline const real3d cone_base_dir_x( void ) const;
        inline const real3d cone_base_dir_y( void ) const { return -direction(); }  // from endpoint() to origin()

        real3d A;
        real3d B;
        real3d B_inv;
        real3d B_norm;
        RAY_KIND KIND;
        real64 SOLID_ANGLE;
        real64 CONE_ANGLE;  
        int top;
    };

    class Polygon
    {
    public:
        uint        mtl_i;                  // index into materials array
        uint        vtx_cnt;                // number of vertices
        uint        vtx_i;                  // index into vertexes array of first vertex
        real3       normal;                 // surface normal
        real        area;                   // surface area
        
        bool bounding_box( const Model * model, AABB& box, real padding=0 ) const;
        bool vertex_info( const Model * model, uint vtx_cnt, real3 p[], real3 n[], real2 uv[] ) const;
        bool vertex_info( const Model * model, uint vtx_cnt, real3 p[] ) const;
        bool sample( const Model * model, HitRecord& rec, real solid_angle=0.0, real dir_sqr=0.0 ) const;  // rec.p must be set on input
        bool hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                  real solid_angle, real t_min, real t_max, HitRecord& rec ) const;
        inline bool is_emissive( const Model * model ) const;
    };

    class Vertex
    {
    public:
        uint        v_i;                    // index into positions array
        uint        vn_i;                   // index into normals array
        uint        vt_i;                   // index into texcoords array
    };

    // see https://www.fileformat.info/format/material/
    // we don't yet support advanced features of the spec
    //
    class Material
    {
    public:
        uint            name_i;             // index into strings array
        real3           Ka;                 // RGB ambient reflectivity                         (default: [1,1,1])
        real3           Kd;                 // RGB diffuse reflectivity                         (default: [1,1,1])
        real3           Ke;                 // RGB emissive                                     (default: [0,0,0])
        real3           Ks;                 // RGB specular reflectivity                        (default: [0,0,0])
        real3           Tf;                 // RGB tranmission filter                           (default: [1,1,1])
        real            Tr;                 // transparency (currently used for RGB)            (default: 0=opaque)
        real            Ns;                 // specular exponent (focus of specular highlight)  (default: 0)
        real            Ni;                 // optical density == index of refraction (0.001 to 10, 1=default=no bending, 1.5=glass)
        real            d;                  // dissolve                                         (default: 1=opaque)
        real            illum;              // ilumination model (0-10, consult documentation)  (default: 2)
        uint            map_Ka_i;           // ambient texture (multiplied by Ka) 
        uint            map_Kd_i;           // diffuse texture (multiplied by Kd)
        uint            map_Ke_i;           // emittance texture (multipled by Ke)
        uint            map_Ks_i;           // specular color texture (multipled by Ks)
        uint            map_Ns_i;           // specular highlight texture (multipled by Ns)
        uint            map_d_i;            // alpha texture (multiplied by d)
        uint            map_Bump_i;         // bump map texture
        uint            map_refl_i;         // reflection map
        uint64_t        unused[8];          // leave room for temporary hacks

        inline bool is_emissive( void ) const { return Ke.c[0] != 0.0 || Ke.c[1] != 0.0 || Ke.c[2] != 0.0; }
    };

    // Graphs will replace simple Materials
    // Work-in-progress:
    //
    enum class GRAPH_NODE_KIND
    {
        BLOCK,                                  // there is currently only one block of assigns
        ASSIGN,                                 // the only statement
        REAL,                                   // one real number constant
        REAL3,                                  // vector constant
        STR,                                    // string constant
        ID,                                     // identifier 
        OP,                                     // operator
        NARY,                                   // n-ary operator or function call (first child is ID or OP, rest are arguments)
    };

    class Graph_Node
    {
    public:
        GRAPH_NODE_KIND kind;
        union {
            uint        child_first_i;          // index of first child node
            uint        s_i;                    // string index (for ID, OP)
            real        r;                      // real constant
            real3       r3;                     // real3 constant
        } u;
        uint            parent_i;
        uint            sibling_i;

        std::string str( const Model * model, std::string indent="" ) const;       // recursive
    };

    uint gph_node_alloc( Model::GRAPH_NODE_KIND kind, uint parent_i=uint(-1) );

    // Volume - made up of one or more VolumeGrid
    //
    // Note: Currently, all grids must have the same bounding boxes.  This is so that
    //       we can return a Volume when a hit() occurs for whatever reason, though
    //       we will also return the grid_i that caused the hit() to stop.
    //

    // VolumeGrid - made up of one or more VolumeGrid
    //
    enum class VolumeVoxelType : uint32_t
    {
        UNKNOWN = 0, 
        FLOAT   = 1, 
        DOUBLE  = 2, 
        INT16   = 3,
        INT32   = 4, 
        INT64   = 5, 
        VEC3F   = 6, 
        VEC3D   = 7,
        MASK    = 8, 
        FP16    = 9
    };

    enum class VolumeVoxelClass : uint32_t
    {
        UNKNOWN      = 0, 
        DENSITY      = 1,                       // particle density
        TEMPERATURE  = 2,                       // used for fire
        RGB          = 3,                       // scattering albedo
        EMISSIVE_RGB = 4,                       // used for lightning
        NORMAL       = 5,                       // surface normal
        IOR          = 6,                       // index of refraction
        F0           = 7,                       // more general 3D version of IOR (see early comments on how to convert from IOR to F0)
        G            = 8,                       // scattering phase function mean cosine (0 == omnidirectional, 1 == forward)
        ATTENUATION  = 9,                       // 3D attenuation (default is 1,1,1 which means none)
        MFPL         = 10                       // mean free path length
    };

    enum class VolumeGridClass : uint32_t
    {
        UNKNOWN      = 0, 
        LEVELSET     = 1, 
        FOGVOLUME    = 2, 
        STAGGERED    = 3,
        POINTINDEX   = 4,
        POINTDATA    = 5
    };

    class Volume
    {
    public:
        uint            name_i;                 // index in strings array (null-terminated strings)
        uint            grid_cnt;               // number of grids
        uint            grid_i;                 // index into volume_grids[] array of first grid

        bool bounding_box( const Model * model, AABB& box, real padding=0 ) const;
        bool hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                  real solid_angle, real t_min, real t_max, HitRecord& rec ) const;
        uint voxel_class_grid_i( const Model * model, VolumeVoxelClass c ) const;
        bool is_emissive( const Model * model ) const;
        bool rand_emissive_xyz( const Model * model, _int& x, _int& y, _int& z ) const;  // returns true and xyz of some emissive voxel, else false if none found

        std::string str( const Model * model, std::string indent="" ) const;       // recursive
    };

    class VolumeGrid
    {
    public:
        uint             name_i;                // index in strings[] array (null-terminated strings)
        uint64           voxel_cnt;             // number of active voxels
        VolumeVoxelType  voxel_type;            // voxel data type (float, double, etc.)
        VolumeVoxelClass voxel_class;           // DENSITY, TEMPERATURE, etc.
        VolumeGridClass  grid_class;            // class of grid (level set, fog, staggered)
        uint             node_cnt[4];           // not sure yet what this is for
        AABB             world_box;             // AABB in world space (what we normally think of as the AABB)
        AABBI            index_box;             // AABB in index space
        real3d           world_voxel_size;      // 3D size of voxel in world units (a voxel is not always a cube)
        real3d           world_voxel_size_inv;  // (1,1,1) / world_voxel_size
        real3d           world_translate;       // 3D translation of scaled xyz index to get to world point
        uint64           voxel_i;               // index into voxels[] of first voxel 

        bool   bounding_box( const Model * model, AABB& box, real padding=0 ) const;
        bool   bounding_box( const Model * model, _int x, _int y, _int z, AABBD& box ) const;
        bool   is_active( const Model * model, _int x, _int y, _int z ) const;
        const void * value_ptr( const Model * model, _int x, _int y, _int z, bool * is_active=nullptr ) const;
        _int   int_value( const Model * model, _int x, _int y, _int z ) const;
        _int64 int64_value( const Model * model, _int x, _int y, _int z ) const;
        real   real_value( const Model * model, _int x, _int y, _int z ) const;
        real64 real64_value( const Model * model, _int x, _int y, _int z ) const;
        real3  real3_value( const Model * model, _int x, _int y, _int z ) const;
        real3d real3d_value( const Model * model, _int x, _int y, _int z ) const;
    };

    class Texture
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        TEXTURE_COMPRESSION compression;    // NONE or ASTC 
        uint            width;              // level 0 width
        uint            height;             // level 0 height
        uint            nchan;              // number of channels (typically 3 or 4)
        uint64          texel_i;            // index into texels array of first texel (uncompressed) or ASTC header
        real            bump_multiplier;    // should be for bump maps only but also used with other textures
        real3           offset;             // uvw offset of texture map on surface (w is 0 for 2D textures)
        real3           scale;              // uvw scale  of texture map on surface (w is 0 for 2D textures)
        uint64_t        unused[8];          // leave room for temporary hacks

        // for ARM .astc files and our internal format
        struct ASTC_Header                      
        {
            uint8_t magic[4];
            uint8_t blockdim_x;
            uint8_t blockdim_y;
            uint8_t blockdim_z;
            uint8_t xsize[3];                       // x-size = xsize[0] + xsize[1] + xsize[2]
            uint8_t ysize[3];                       // x-size, y-size and z-size are given in texels;
            uint8_t zsize[3];                       // block count is inferred

            // utility functions
            bool        magic_is_good( void ) const;
            uint        size_x( void ) const;                 // width of texture etc.
            uint        size_y( void ) const;
            uint        size_z( void ) const;
            uint        blk_cnt_x( void ) const;              
            uint        blk_cnt_y( void ) const;
            uint        blk_cnt_z( void ) const;
            uint        blk_cnt( void ) const;                // including header block
            uint        byte_cnt( void ) const;               // blk_cnt() * 16
        };

        real4           texel_read(              const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr=nullptr, uint64 * byte_cnt=nullptr,
                                                 bool do_srgb_to_linear=false ) const;
        static void     astc_blk_dims_get( uint width, uint height, uint& blk_dim_x, uint& blk_dim_y );

        // self-contained static function for implementing mip-based addressing 
        static void     astc_blk_addr_get( uint      ray_id,
                                           uint      blockdim_x, uint blockdim_y, uint blockdim_z,
                                           uint      xsize,      uint ysize,      uint zsize,
                                           real      u,          real v,          real w,
                                           real      frac_uv_covg, 
                                           uint      size_w,
                                           uint      addr_w,
                                           uint      uv_frac_w,
                                           uint      covg_frac_w,
                                           uint64_t  blks_addr,
                                           uint64_t& blk_addr, uint& s, uint& t, uint& r );

        // self-contained static function for decoding a single texel in a block whose data is supplied
        // note: only blockdim* need to be set in the ASTC_Header (nothing else is used or checked)
        static real4    astc_decode_texel( const unsigned char * bdata, const ASTC_Header * astc_hdr, 
                                           uint s, uint t, uint r, bool do_srgb_to_linear );

    private:
        real4           texel_read_uncompressed( const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const;
        real4           texel_read_astc(         const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const;

        static void     astc_decode_block_mode( const unsigned char * bdata, unsigned char * rbdata, uint& plane_cnt, uint& partition_cnt, 
                                                uint& weights_w, uint& weights_h, uint& weights_d, uint& R, uint& H, uint& weights_qmode );
        static uint     astc_decode_partition( const unsigned char * bdata, uint x, uint y, uint z, uint partition_cnt, uint texel_cnt );
        static void     astc_decode_color_endpoint_mode( const unsigned char * bdata, uint plane_cnt, uint partition, uint partition_cnt, 
                                                         uint weights_bit_cnt, uint below_weights_start,
                                                         uint& cem, uint& qmode, uint& trits, uint& quints, uint& bits, 
                                                         uint& endpoints_cnt, uint& endpoints_start, uint& endpoints_v_first, uint& endpoints_v_cnt );
        static void     astc_decode_color_endpoints( const unsigned char * bdata, 
                                                     uint cem, uint cem_qmode, uint trits, uint quints, uint bits, 
                                                     uint endpoints_cnt, uint endpoints_start, uint endpoints_v_first, uint endpoints_v_cnt,
                                                     uint endpoints[4][2] );
        static void     astc_decode_bit_transfer_signed( int& a, int& b );
        static void     astc_decode_blue_contract( int& r, int& g, int& b, int& a );
        static uint     astc_decode_clamp_unorm8( int v );
        static void     astc_decode_weights( const unsigned char * rbdata, uint weights_start, uint plane_cnt,
                                             uint blk_w, uint blk_h, uint blk_d, uint weights_w, uint weights_h, uint weights_d,
                                             uint s, uint t, uint r, uint trits, uint quints, uint bits, uint qmode,
                                             uint plane_weights[2] );
        static uint     astc_decode_integer( const unsigned char * bdata, uint start, uint cnt, uint i, uint trits, uint quints, uint bits );
        static std::string astc_dat_str( const unsigned char * bdata ); 


        struct ASTC_Range_Encoding
        {
            uint max_p1;        // max value + 1
            uint trits;         // if 1, MSBs are trit-encoded
            uint quints;        // if 1, MSBs are quint-encoded
            uint bits;          // number of unshared LSBs for value
        };

        static const ASTC_Range_Encoding astc_range_encodings[21];        
        static const uint                astc_color_quantization_mode[17][128];   // [color_int_cnt/2][color_bits]

        static const uint                astc_integer_trits[256][5];
        static const uint                astc_integer_quints[128][3];

        static const uint                astc_weight_unquantized[12][32];         // [quantization_mode][qv]
        static const uint                astc_color_unquantized[21][256];         // [quantization_mode][qv]

        struct ASTC_Texel_Decimation_Encoding  // per texel in a block
        {
            uint        weight_cnt;         // number of weights to use; up to 4 max
            uint        weight_i[4];        // index of weight in weights[] array
            uint        weight_factor[4];   // factor to multiply by
        };
        using ASTC_Block_Weights_Decimation_Encoding = std::vector<ASTC_Texel_Decimation_Encoding>;         // per block dims by weight dims
        using ASTC_Block_Decimation_Encoding         = std::vector<ASTC_Block_Weights_Decimation_Encoding>; // per block dims
        using ASTC_Decimation_Encoding               = std::vector<ASTC_Block_Decimation_Encoding>;         // for all possible blocks
        static ASTC_Decimation_Encoding   astc_weight_decimation_encodings;                                 // this gets filled in on demand by...

        // the following will fill in the table based on demand for certain combinations of block and weights dims
        static const ASTC_Block_Weights_Decimation_Encoding& astc_block_weights_decimation_encoding( uint blk_w, uint blk_h, uint weights_w, uint weight_h );  
    };

    enum class BVH_NODE_KIND
    {
        BVH_NODE,                           // BVH_Node
        POLYGON,                            // Polygon
        INSTANCE,                           // Instance
        VOLUME,                             // Volume
    };

    class BVH_Node
    {
    public:
        AABB            box;                // bounding box
        BVH_NODE_KIND   left_kind;          // see kinds above
        BVH_NODE_KIND   right_kind;         // see kinds above
        uint            left_i;             // index into appropriate kind array for left subtree 
        uint            right_i;            // index into appropriate kind array for right subtree 

        bool bounding_box( const Model * model, AABB& b ) const;
        bool hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                  real solid_angle, real t_min, real t_max, HitRecord& rec ) const;
    };

    class Matrix                            // 4x4 transformation matrix used by Camera and Instance
    {
    public:
        real            m[4][4];

        inline Matrix( void )          { identity(); }
        static Matrix make_view( const real3& lookfrom, const real3& lookat, const real3& vup );
        static Matrix make_frustum( real left, real right, real bottom, real top, real near, real far );
        static Matrix make_perspective( real vfov, real aspect, real near, real far );

        bool   is_identity( void ) const;               // returns true if this is the identity matrix (can speed up transform, etc.)
        void   identity(  void );                       // make this the identity matrix
        void   translate( const real3& translation );   // translate this matrix by a real3
        void   scale(     const real3& scaling );       // scale this matrix by a real3
        void   rotate_xz( double radians );             // rotate by radians in xz plane (yaw)
        void   rotate_yz( double radians );             // rotate by radians in yz plane (pitch)
        void   rotate_xy( double radians );             // rotate by radians in xy plane (roll)

        Matrix operator + ( const Matrix& m ) const;    // add two matrices
        Matrix operator - ( const Matrix& m ) const;    // subtract two matrices
        bool   operator == ( const Matrix& m ) const;   // return true if both matrices are equal
        void   multiply( double s );                    // multiply this matrix by scalar

        real4  row( uint r ) const;                     // returns row r as a vector
        real4  column( uint c ) const;                  // returns column c as a vector
        void   transform( const real4& v, real4& r ) const; // r = *this * v
        void   transform( const real3& v, real3& r, bool div_by_w=true ) const;   // r = *this * v (and divide by w by default)
        void   transform( const real3d& v, real3d& r, bool div_by_w=true ) const; 
        void   transform( const Ray& r, Ray& r2 ) const;
        void   transform( const Ray64& r, Ray64& r2 ) const;
        void   transform( const Matrix& M2, Matrix& M3 ) const; // M3 = *this * M2
        void   transpose( Matrix& mt ) const;           // return the transpose this matrix 
        void   invert( Matrix& minv ) const;            // return the inversion this matrix
        void   invert_affine( Matrix& minv ) const;     // same but faster because assumes an affine transform
        void   adjoint( Matrix& M ) const;              // used by invert()
        double determinant( void ) const;               // returns the determinant (as double for high-precision) 
        void   cofactor( Matrix& C ) const;             // used by adjoint() and determinant() 
        double subdeterminant( uint exclude_row, uint exclude_col ) const; // used by cofactor() 
    };

    class Camera 
    {
    public:
        Camera( real3 lookfrom, real3 lookat, real3 vup, real vfov, real near, real far, real aspect, real aperture, real focus_dist, int nx, int ny, int spp, uint name_i=uint(-1) );

        Ray get_ray( real s_pixel, real t_pixel, real s_lens, real t_lens ) const;  // simple thin lens camera by default
        bool overlaps_box( const AABB& box ) const;                                 // camera-box intersection test

        // inputs
        real3           lookfrom;           // camera location
        real3           lookat;             // point camera is aimed at
        real3           vup;                // camera up direction
        real            aperture;           // lens aperture
        real            focus_dist;         // aka focal_length
        real            near;               // near clip plane distance (typical: 0.1)
        real            far;                // far clip plane distance (typical: 10000)
        real            vfov;               // vertical field of view
        real            aspect;             // aspect ratio
        uint            name_i;             // index in strings array (null-terminated strings)

        // derived
        real3           lower_left_corner;
        real3           horizontal;
        real3           vertical;
        real3           u, v, w;
        real32          lens_radius;
        real32          solid_angle;
        Matrix          perspective_projection;

        struct Plane {
            real3       normal;             // [a, b, c] from standard unnormalized plane equation 
            real        distance;           // d         from standard unnormalized plane equation
        };
        Plane           frustum_planes[6];  // 0=near 1=far 2=bottom 3=top left=4 right=5

    private:
        void concentric_point_on_unit_disk(real s, real t, real& x, real& y) const;
    };

    enum class INSTANCE_KIND
    {
        MODEL,                              // instance is another Model (by name)
        MODEL_PTR,                          // instance is another Model (read in and with a resolved pointer)
        OBJECT,                             // instance is an Object within the current model
    };

    class Instance      
    {
    public:
        INSTANCE_KIND   kind;               // see above
        uint            name_i;             // index in strings array of instance name
        uint            model_name_i;       // index in strings array of name of instanced Model 
        uint            model_file_name_i;  // index in strings array of file name of instanced Model 
        uint            matrix_i;           // index of Matrix in matrixes[] array (transformation matrix)
        uint            matrix_inv_i;       // index of inverted Matrix in matrixes[] array (inverse transformation matrix)
        uint            matrix_inv_trans_i; // index of transposed inverted Matrix in matrixes[] array (used for transforming normals)
        AABB            box;                // bounding box in global world space, not instance space
        union {
            Model *     model_ptr;          // resolved Model pointer for MODEL_PTR
            uint        obj_i;              // index of Object in objects[] array
        } u;

        bool bounding_box( const Model * model, AABB& box, real padding=0 ) const;
        bool hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                  real solid_angle, real t_min, real t_max, HitRecord& rec );
    };

    enum class LIGHT_KIND
    {
        DIRECTIONAL,
        POINT,
        SPHERICAL,
        RECTANGULAR,
        PROBE,
    };

    class Light
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        uint            file_name_i;        // index in strings array of file name for light probe only
        LIGHT_KIND      kind;               // type of light
        real3           intensity;          // aka color
        real3           direction;          // point and directional lights only
        real3           position;           // point and spherical lights only
        real3           extent;             // spherical (r, r, r) and rectangular (w, h, 0) lights only
        real            opening_angle;      // point lights only
        real            penumbra_angle;     // point lights only
        uint            diffuse_sample_cnt; // applies to light probe only
        uint            specular_sample_cnt;// ditto
        uint64_t        unused[8];          // leave room for temporary hacks
    };

    class Frame                             // one frame in an animation sequence
    {
    public:
        real            time;               // in seconds (divide by Header.animation_speed)
        uint            camera_i;           // index into cameras array (name = animation name + frame index)
    };

    class Animation
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        uint            start_camera_i;     // index in cameras array of initial camera to use
        uint            start_frame_i;      // index in animations array of first frame in animation sequence
        uint            frame_cnt;          // number of frames in animation sequence
        bool            is_looped;          // whether the animation sequence should be looped (default: false)
    };

    // useful methods for allocating above structures in this Model's arrays with proper resizing, etc.
    // these all return indexes, which is what the user should store as pointers can change as arrays are resized later
    inline uint make_string( std::string s );     // returns index in strings[] array
    inline uint make_texture( std::string name, std::string dir );                                                      // from file 
    inline uint make_texture( std::string name, const unsigned char * data, uint w, uint h, uint nchan );               // from w * h * nchan bytes already in memory
    inline uint make_constant_texture( std::string name, const real4& value );                                          // from one RGBA value 
    inline uint make_constant_texture( std::string name, const real3& value );                                          // from one RGB  value 
    inline uint make_material( std::string name, real3 Ka, real3 Kd, real3 Ke, real3 Ks, real3 Tf, real Tr, real Ns, real Ni, real d, real illum,
                               uint map_Ka_i, uint map_Kd_i, uint map_Ke_i, uint map_Ks_i, uint map_Ns_i, uint map_d_i, uint map_Bump_i, uint map_refl_i );
    inline uint make_emissive_material( std::string name, real3 Ke, uint map_Ke_i=uint(-1), real d=1.0 );                        
    inline uint make_diffuse_material( std::string name, real3 Kd, uint map_Kd_i=uint(-1), real d=1.0, uint map_d_i=uint(-1) );

    inline uint make_vertex( const real3& p, const real3& n, const real2& uv );                                         // position, normal, and texcoord
    inline uint make_polygon( uint vtx_cnt, const real3 p[], const real3 n[], const real2 uv[], uint mtl_i=uint(-1) );  // positions, normals, texcoords
    inline uint make_polygon( uint vtx_cnt, const real3 p[], uint mtl_i=uint(-1) );                                     // position, use surface normal, fudge texcoords
    inline uint make_sphere( const real3& center, real radius, uint nlong, uint nlat, uint mtl_i=uint(-1) );            // center, radius, long+lat divisions
    inline uint make_box( const real3 front[4], const real3 back[4], uint mtl_i=uint(-1) );                             // front and back face vertices (CCW from bottom-left as looking from front)
    inline uint make_aabox( const AABB& box, uint mtl_i=uint(-1) );                                                     // axis-aligned box (simpler case)
    inline uint make_skybox( const real3& center, real radius, const real3& up, uint nlong, uint nlat, uint mtl_i=uint(-1) ); // center, radius, up vector, long+lat divisions
    inline uint make_matrix( void );                                                                                    // returns an identity matrix for starters
    inline uint make_instance( std::string name, std::string model_name, Model * model, uint matrix=uint(-1), std::string file_name="" ); // MODEL_PTR instance 

    // call this to rebuild the BVH after adding all primitives dynamically
    void                bvh_build( BVH_TREE bvh_tree=BVH_TREE::BINARY );

    // srgb<->linear conversion utilities
    static real srgb_to_linear_gamma( real srgb );
    static real srgb8_to_linear_gamma( uint8_t srgb );   // more efficient, uses LUT
    static real linear_to_srgb_gamma( real linear );
    static real linear8_to_srgb_gamma( uint8_t linear ); // more efficient, uses LUT

    // system utilities
    void cmd( std::string c, std::string error="command failed", bool echo=true );  // calls std::system and aborts if not success

    // memory allocation utilities
    template<typename T> inline T *  aligned_alloc( uint64 cnt );
    template<typename T> inline void perhaps_realloc( T *& array, const uint64& hdr_cnt, uint64& max_cnt, uint64 add_cnt );

    // file utilities
    static void dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); 
    static bool file_exists( std::string file_name );
    bool file_read( std::string file_name, char *& start, char *& end, bool is_read_only=true );  // sucks in entire file
    bool file_read( std::string file_name, const char *& start, const char *& end );    
    void file_write( std::string file_name, const unsigned char * data, uint64_t byte_cnt );

    // parsing utilities for files sucked into memory
    uint line_num;
    bool skip_whitespace_to_eol( const char *& xxx, const char * xxx_end );  // on this line only
    bool skip_whitespace( const char *& xxx, const char * xxx_end );
    bool skip_to_eol( const char *& xxx, const char * xxx_end );
    bool skip_through_eol( const char *& xxx, const char * xxx_end );
    bool eol( const char *& xxx, const char * xxx_end );
    bool expect_char( char ch, const char *& xxx, const char * xxx_end, bool skip_whitespace_first=false );
    bool expect_eol( const char *& xxx, const char * xxx_end );
    bool expect_cmd( const char * s, const char *& xxx, const char * xxx_end );
    bool parse_string( std::string& s, const char *& xxx, const char * xxx_end );
    bool parse_string_i( uint& s, const char *& xxx, const char * xxx_end );
    bool parse_name( const char *& name, const char *& xxx, const char * xxx_end );
    bool parse_id( std::string& id, const char *& xxx, const char * xxx_end );
    bool parse_id_i( uint& id_i, const char *& xxx, const char * xxx_end );
    bool parse_option_name( std::string& option_name, const char *& xxx, const char * xxx_end );
    bool parse_real3( real3& r3, const char *& xxx, const char * xxx_end, bool has_brackets=false );
    bool parse_real2( real2& r2, const char *& xxx, const char * xxx_end, bool has_brackets=false );
    bool parse_real( real& r, const char *& xxx, const char * xxx_end, bool skip_whitespace_first=false );
    bool parse_real64( real64& r, const char *& xxx, const char * xxx_end, bool skip_whitespace_first=false );
    bool parse_int( _int& i, const char *& xxx, const char * xxx_end );
    bool parse_int64( _int64& i, const char *& xxx, const char * xxx_end );
    bool parse_uint( uint& u, const char *& xxx, const char * xxx_end, uint base=10 );
    bool parse_uint64( uint64& u, const char *& xxx, const char * xxx_end, uint base=10 );
    bool parse_bool( bool& b, const char *& xxx, const char * xxx_end );
    std::string surrounding_lines( const char *& xxx, const char * xxx_end );

    // structs
    void *              mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    size_t              mapped_region_len;  
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays
    char *              strings;
    Object *            objects;
    Polygon *           polygons;
    uint *              emissive_polygons;  
    Vertex *            vertexes;
    real3 *             positions;
    real3 *             normals;
    real2 *             texcoords;
    Material *          materials;
    Graph_Node *        graph_nodes;
    Volume *            volumes;
    VolumeGrid *        volume_grids;
    unsigned char *     voxels;
    Texture *           textures;
    unsigned char *     texels;
    BVH_Node *          bvh_nodes;
    Matrix *            matrixes;
    Instance *          instances;
    Light *             lights;
    Camera *            cameras;
    Frame *             frames;
    Animation *         animations;

    // maps of names to array indexes
    std::map<std::string, uint>    name_to_obj_i;
    std::map<std::string, uint>    name_to_mtl_i;
    std::map<std::string, uint>    name_to_tex_i;
    std::map<std::string, uint>    name_to_light_i;
    std::map<std::string, Model *> name_to_model;               // file path name to instanced Model
    std::map<std::string, uint>    name_to_camera_i;
    std::map<std::string, uint>    name_to_animation_i;

    // map of model file name to ignored polygon indexes
    std::map<std::string, std::vector<uint> *> model_ignored_polys;

    // debug flags
    static bool         debug;
    static uint         debug_tex_i;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//
//  ___                 _                           _        _   _
// |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
//  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
//  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
// |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
//               |_|
//
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
private:
    typedef enum 
    {
        CMD_O,
        CMD_G,
        CMD_S,
        CMD_V,
        CMD_VN,
        CMD_VT,
        CMD_F,
        CMD_MTLLIB,
        CMD_USEMTL
    } obj_cmd_t;
            
    typedef enum 
    {
        CMD_NEWMTL,
        CMD_KA,
        CMD_KD,
        CMD_KE,
        CMD_KS,
        CMD_TF,
        CMD_TR,
        CMD_NS,
        CMD_NI,
        CMD_D,
        CMD_ILLUM,
        CMD_MAP_KA,
        CMD_MAP_KD,
        CMD_MAP_KE,
        CMD_MAP_KS,
        CMD_MAP_NS,
        CMD_MAP_D,
        CMD_MAP_BUMP,
        CMD_MAP_REFL
    } mtl_cmd_t;
            
    const char * fsc_start;
    const char * fsc_end;
    const char * fsc;
    const char * obj_start;
    const char * obj_end;
    const char * obj;
    const char * mtl_start;
    const char * mtl_end;
    const char * mtl;
    const char * gph_start;
    const char * gph_end;
    const char * gph;
    const char * nvdb_start;
    const char * nvdb_end;
    const char * nvdb;

    const uint64_t **   volume_to_emissive;             // information about emissive voxels

    bool load_fsc( std::string fsc_file, std::string dir_name );        // .fscene 
    bool load_obj( std::string obj_file, std::string dir_name );        // .obj
    bool load_mtl( std::string mtl_file, std::string dir_name, bool replacing=false );
    bool load_tex( const char * tex_name, std::string dir_name, Texture *& texture, const unsigned char * data=nullptr, uint w=0, uint h=0, uint nchan=0 );
    bool load_gph( std::string gph_file, std::string dir_name, std::string base_name, std::string ext_name );
    bool load_nvdb( std::string nvdb_file, std::string dir_name, std::string base_name, std::string ext_name );

    bool parse_obj_cmd( obj_cmd_t& cmd );
    bool parse_mtl_cmd( mtl_cmd_t& cmd, const char *& mtl, const char * mtl_end );
    bool parse_gph_real( uint& r_i, const char *& gph, const char * gph_end );
    bool parse_gph_real3( uint& r3_i, const char *& gph, const char * gph_end );
    bool parse_gph_str( uint& s_i, const char *& gph, const char * gph_end );
    bool parse_gph_id( uint& id_i, const char *& gph, const char * gph_end );
    bool parse_gph_op_lookahead( uint& s_i, const char *& gph, const char * gph_end );
    bool parse_gph_op( uint& op_i, const char *& gph, const char * gph_end );
    bool parse_gph_term( uint& term_i, const char *& gph, const char * gph_end );
    bool parse_gph_expr( uint& expr_i, const char *& gph, const char * gph_end, uint curr_prec=0 );

    // BVH builder
    uint bvh_qsplit( BVH_NODE_KIND kind, uint poly_i, uint n, real pivot, uint axis );
    uint bvh_node( BVH_NODE_KIND kind, uint i, uint n, uint axis );

    bool write_uncompressed( std::string file_path );
    bool read_uncompressed( std::string file_path );
};

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// USEFUL CONSTANTS (IN FLOAT)
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
constexpr Model::real32 PI          = M_PI;
constexpr Model::real32 PI2         = M_PI*2.0;          // aka TAU
constexpr Model::real32 PI4         = M_PI*4.0;
constexpr Model::real32 ONE_OVER_PI = 1.0/M_PI;
constexpr Model::real32 ONE_OVER_PI2= 1.0/(2.0*M_PI);    
constexpr Model::real32 ONE_OVER_PI4= 1.0/(4.0*M_PI);    
constexpr Model::real32 PI_DIV_2    = M_PI/2.0;
constexpr Model::real32 PI_DIV_4    = M_PI/4.0;
constexpr Model::real32 PI_DIV_180  = M_PI/180.0;
constexpr Model::real32 PI_DIV_360  = M_PI/360.0;

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// USEFUL REAL FUNCTIONS
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real   mulby2( const Model::real& x )                                  { return 2.0 * x;                }
inline Model::real64 mulby2( const Model::real64& x )                                { return 2.0 * x;                }
inline Model::real   mulby4( const Model::real& x )                                  { return 4.0 * x;                }
inline Model::real64 mulby4( const Model::real64& x )                                { return 4.0 * x;                }
inline Model::real   divby2( const Model::real& x )                                  { return 0.5 * x;                }
inline Model::real64 divby2( const Model::real64& x )                                { return 0.5 * x;                }
inline Model::real   divby4( const Model::real& x )                                  { return 0.25 * x;               }
inline Model::real64 divby4( const Model::real64& x )                                { return 0.25 * x;               }
inline Model::real   sqr( const Model::real& x )                                     { return x*x;                    }
inline Model::real64 sqr( const Model::real64& x )                                   { return x*x;                    }
inline Model::real   rcp( const Model::real& x )                                     { return 1.0 / x;                }
inline Model::real64 rcp( const Model::real64& x )                                   { return 1.0 / x;                }
inline Model::real   rsqrt( const Model::real& x )                                   { return 1.0 / std::sqrt( x );   }
inline Model::real64 rsqrt( const Model::real64& x )                                 { return 1.0 / std::sqrt( x );   }
inline Model::real   clamp( const Model::real& x )                                   { return std::fmin(std::fmax(x, 0.0), 1.0);} 
inline Model::real64 clamp( const Model::real64& x )                                 { return std::fmin(std::fmax(x, 0.0), 1.0);} 
inline Model::real   clamp( const Model::real& x, const Model::real& min, const Model::real& max )   { return std::fmin(std::fmax(x, min), max);} 
inline Model::real64 clamp( const Model::real64& x, const Model::real64& min, const Model::real64& max )   { return std::fmin(std::fmax(x, min), max);} 
inline Model::real   lerp( const Model::real& x, const Model::real& y, const Model::real& t )        { return (1.0 - t)*x + t*y;      }
inline Model::real64 lerp( const Model::real64& x, const Model::real64& y, const Model::real64& t )  { return (1.0 - t)*x + t*y;      }
inline Model::real   expc( double b, const Model::real& x )                          { return std::pow( b, x );       }
inline Model::real64 expc( double b, const Model::real64& x )                        { return std::pow( b, x );       }
inline Model::real   hypoth( const Model::real& x, const Model::real& y )            { return std::sqrt( x*x - y*y ); }
inline Model::real64 hypoth( const Model::real64& x, const Model::real64& y )        { return std::sqrt( x*x - y*y ); }
inline Model::real   hypoth1( const Model::real& y )                                 { return std::sqrt( 1.0 - y*y ); }
inline Model::real64 hypoth1( const Model::real64& y )                               { return std::sqrt( 1.0 - y*y ); }
inline void sincos( const Model::real& x, Model::real& si, Model::real& co )         { si = std::sin( x ); co = std::cos( x ); }
inline void sincos( const Model::real& x, Model::real& si, Model::real& co, const Model::real& r ) { si = r*std::sin( x ); co = r*std::cos( x ); }
inline void sincos( const Model::real64& x, Model::real64& si, Model::real64& co )   { si = std::sin( x ); co = std::cos( x ); }
inline void sincos( const Model::real64& x, Model::real64& si, Model::real64& co, const Model::real64& r ) { si = r*std::sin( x ); co = r*std::cos( x ); }

inline static void minmax( const Model::real x0, const Model::real x1, const Model::real x2, Model::real &min, Model::real &max )
{
    min = std::min(std::min(x0, x1), x2);
    max = std::max(std::max(x0, x1), x2);
}
inline static void minmax( const Model::real64 x0, const Model::real64 x1, const Model::real64 x2, Model::real64 &min, Model::real64 &max )
{
    min = std::min(std::min(x0, x1), x2);
    max = std::max(std::max(x0, x1), x2);
}

inline Model::real length_2D_line(Model::real x0, Model::real y0, Model::real x1, Model::real y1) {
    return std::sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
}
inline Model::real64 length_2D_line(Model::real64 x0, Model::real64 y0, Model::real64 x1, Model::real64 y1) {
    return std::sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
}

inline Model::real area_2D_triangle(Model::real x0, Model::real y0, Model::real x1, Model::real y1, Model::real x2, Model::real y2) {
    Model::real A = length_2D_line(x0,y0,x1,y1);
    Model::real B = length_2D_line(x0,y0,x2,y2);
    Model::real C = length_2D_line(x2,y2,x1,y1);
    Model::real S = 0.5*(A+B+C);
    return std::sqrt(S*(S-A)*(S-B)*(S-C));
}
inline Model::real64 area_2D_triangle(Model::real64 x0, Model::real64 y0, Model::real64 x1, Model::real64 y1, Model::real64 x2, Model::real64 y2) {
    Model::real64 A = length_2D_line(x0,y0,x1,y1);
    Model::real64 B = length_2D_line(x0,y0,x2,y2);
    Model::real64 C = length_2D_line(x2,y2,x1,y1);
    Model::real64 S = 0.5*(A+B+C);
    return std::sqrt(S*(S-A)*(S-B)*(S-C));
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// USEFUL VECTOR FUNCTIONS
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real2  operator +( const Model::real2& v1, const Model::real2& v2 )                               { return Model::real2(v1.c[0] + v2.c[0], v1.c[1] + v2.c[1]);            }
inline Model::real2  operator -( const Model::real2& v1, const Model::real2& v2 )                               { return Model::real2(v1.c[0] - v2.c[0], v1.c[1] - v2.c[1]);            }
inline Model::real2  operator *( const Model::real2& v1, const Model::real2& v2 )                               { return Model::real2(v1.c[0] * v2.c[0], v1.c[1] * v2.c[1]);            }
inline Model::real2  operator /( const Model::real2& v1, const Model::real2& v2 )                               { return Model::real2(v1.c[0] / v2.c[0], v1.c[1] / v2.c[1]);            }
inline Model::real2  operator *( Model::real t, const Model::real2& v )                                         { return Model::real2(t*v.c[0], t*v.c[1]);                              }
inline Model::real2  operator /( Model::real2 v, Model::real t )                                                { return Model::real2(v.c[0]/t, v.c[1]/t);                              }
inline Model::real2  operator *( const Model::real2& v, Model::real t )                                         { return Model::real2(t*v.c[0], t*v.c[1]);                              }

inline Model::real2  min( const Model::real2& v0, const Model::real2& v1 )                                      { return Model::real2(std::min(v0[0], v1[0]), std::min(v0[1], v1[1]));  }
inline Model::real2  max( const Model::real2& v0, const Model::real2& v1 )                                      { return Model::real2(std::max(v0[0], v1[0]), std::max(v0[1], v1[1]));  }
inline Model::real2  sqr( const Model::real2& v )                                                               { return Model::real2(::sqr(v.c[0]), ::sqr(v.c[1]));                    }
inline Model::real   dot( const Model::real2& v0, const Model::real2& v1 )                                      { return v0.dot( v1 );                                                  }
inline Model::real2  lerp( const Model::real2& v0, const Model::real2& v1, Model::real t )                      { return (Model::real(1.0) - t) * v0 + t * v1;                          }
inline Model::real2  lerp( const Model::real2 &v0, const Model::real2 &v1, const Model::real2 &vt)              { return (Model::real2(1,1) - vt) * v0 + vt * v1;                       }

inline Model::real3  operator +( const Model::real3& v1, const Model::real3& v2 )                               { return Model::real3(v1.c[0] + v2.c[0], v1.c[1] + v2.c[1], v1.c[2] + v2.c[2]); }
inline Model::real3  operator -( const Model::real3& v1, const Model::real3& v2 )                               { return Model::real3(v1.c[0] - v2.c[0], v1.c[1] - v2.c[1], v1.c[2] - v2.c[2]); }
inline Model::real3  operator *( const Model::real3& v1, const Model::real3& v2 )                               { return Model::real3(v1.c[0] * v2.c[0], v1.c[1] * v2.c[1], v1.c[2] * v2.c[2]); }
inline Model::real3  operator /( const Model::real3& v1, const Model::real3& v2 )                               { return Model::real3(v1.c[0] / v2.c[0], v1.c[1] / v2.c[1], v1.c[2] / v2.c[2]); }
inline Model::real3  operator *( Model::real t, const Model::real3& v )                                         { return Model::real3(t*v.c[0], t*v.c[1], t*v.c[2]);                    }
inline Model::real3  operator /( Model::real3 v, Model::real t )                                                { return Model::real3(v.c[0]/t, v.c[1]/t, v.c[2]/t);                    }
inline Model::real3  operator *( const Model::real3& v, Model::real t )                                         { return Model::real3(t*v.c[0], t*v.c[1], t*v.c[2]);                    }

inline Model::real3  min( const Model::real3& v0, const Model::real3& v1 )                                      { return Model::real3(std::min(v0[0], v1[0]), std::min(v0[1], v1[1]), std::min(v0[2], v1[2])); }
inline Model::real3  max( const Model::real3& v0, const Model::real3& v1 )                                      { return Model::real3(std::max(v0[0], v1[0]), std::max(v0[1], v1[1]), std::max(v0[2], v1[2])); }
inline Model::real3  sqr( const Model::real3& v )                                                               { return Model::real3(::sqr(v.c[0]), ::sqr(v.c[1]), ::sqr(v.c[2]));     }
inline Model::real   dot( const Model::real3& v0, const Model::real3& v1 )                                      { return v0.dot( v1 );                                                  }
inline Model::real3  cross( const Model::real3& v0, const Model::real3& v1 )                                    { return v0.cross( v1 );                                                }
inline Model::real3  lerp( const Model::real3& v0, const Model::real3& v1, Model::real t )                      { return (Model::real(1.0) - t) * v0 + t * v1;                          }
inline Model::real3  lerp( const Model::real3 &v0, const Model::real3 &v1, const Model::real3 &vt)              { return (Model::real3(1,1,1) - vt) * v0 + vt * v1;                     }

inline Model::real3 spherical_to_cartesian( const Model::real3& center, Model::real radius, Model::real theta, Model::real phi ) 
{
    Model::real x = std::cos(phi)*std::sin(theta);
    Model::real y = std::sin(phi)*std::sin(theta);
    Model::real z = std::cos(theta);
    return center + radius*Model::real3(x,y,z);
}

// find texture coordinates for a point p on a unit sphere centered at origin
// u and v both shifted to have minimum at 0 and then scaled to 0,1
// u replacerd with 1-u to ensure proper glones from standard maps
void get_sphere_uv( const Model::real3& p, Model::real& u, Model::real& v) {
    Model::real phi = std::atan2(p.z(), p.x());
    Model::real y = p.y();
    if ( y >  1.0 ) y =  1.0;
    if ( y < -1.0 ) y = -1.0;
    Model::real theta = std::asin(y);
    u = 1.0 - (phi + PI) * ONE_OVER_PI2;
    v = (theta + PI_DIV_2) * ONE_OVER_PI;
}

// more complex than polar, but according to 
// A Realistic Camera Model for Computer Graphics (Kolb et al)
// has 15% improvement in error over polar for camera lens sampling
// possibly similar improvement for sampling cones?
inline void concentric_point_on_unit_disk(Model::real s, Model::real t, Model::real& x, Model::real& y) {
    Model::real a = mulby2(s) - 1.0;
    Model::real b = mulby2(t) - 1.0;
    Model::real rad, theta;

    if (std::abs(a) > std::abs(b)) {
        rad = a;
        theta = PI_DIV_4*(b/a);
    }
    else {
        rad = b;
        theta = PI_DIV_2 - PI_DIV_4*(a/b);
    }
    sincos(theta, y, x, rad);
}

// mapping from square to disk then to cartesian
// Shirley/Chiu JGT concentric mapping
inline void concentric_square_disk(Model::real x, Model::real y, Model::real& u, Model::real& v) {
    Model::real  a = 2.0*x-1.0;
    Model::real  b = 2.0*y-1.0;
    Model::real r, phi;
    if (a > -b) {
        if (a > b) {
            r = a;
            phi = (PI_DIV_4)*(b/a) ;
        }
        else {
            r = b;
            phi = (PI_DIV_4)*(2.0-(a/b));
        }
    }
    else {
        if (a < b) {
            r = -a;
            phi = (PI_DIV_4)*(4.0 + (b/a));
        }
        else {
            r = -b;
            if (b != 0.0)
                phi = (PI_DIV_4)*(6.0 - (a/b));
            else
                phi = 0.0;
        }
    }
    u = r*std::cos(phi);
    v = r*std::sin(phi);
}

// Pixar's JGT paper.   Produces two tangent vectors that
// are an orthonormal basis b1, b2, n
inline void branchlessONB(const Model::real3 &n, Model::real3 &b1, Model::real3 &b2)
{
    Model::real sign = std::copysign(Model::real(1), n[2]);
    const Model::real a = -rcp(sign + n[2]);
    const Model::real b = n[0] * n[1] * a;
    b1 = Model::real3(1.0 + sign * sqr(n[0]) * a, sign * b, -sign * n[0]);
    b2 = Model::real3(b, sign + sqr(n[1]) * a, -n[1]);
}
        
inline Model::real3d operator +( const Model::real3d& v1, const Model::real3d& v2 )                             { return Model::real3d(v1.c[0] + v2.c[0], v1.c[1] + v2.c[1], v1.c[2] + v2.c[2]); }
inline Model::real3d operator -( const Model::real3d& v1, const Model::real3d& v2 )                             { return Model::real3d(v1.c[0] - v2.c[0], v1.c[1] - v2.c[1], v1.c[2] - v2.c[2]); }
inline Model::real3d operator *( const Model::real3d& v1, const Model::real3d& v2 )                             { return Model::real3d(v1.c[0] * v2.c[0], v1.c[1] * v2.c[1], v1.c[2] * v2.c[2]); }
inline Model::real3d operator /( const Model::real3d& v1, const Model::real3d& v2 )                             { return Model::real3d(v1.c[0] / v2.c[0], v1.c[1] / v2.c[1], v1.c[2] / v2.c[2]); }
inline Model::real3d operator *( Model::real64 t, const Model::real3d& v )                                      { return Model::real3d(t*v.c[0], t*v.c[1], t*v.c[2]);                   }
inline Model::real3d operator /( Model::real3d v, Model::real64 t )                                             { return Model::real3d(v.c[0]/t, v.c[1]/t, v.c[2]/t);                   }
inline Model::real3d operator *( const Model::real3d& v, Model::real64 t )                                      { return Model::real3d(t*v.c[0], t*v.c[1], t*v.c[2]);                   }

inline Model::real3d min( const Model::real3d& v0, const Model::real3d& v1 )                                    { return Model::real3d(std::min(v0[0], v1[0]), std::min(v0[1], v1[1]), std::min(v0[2], v1[2])); }
inline Model::real3d max( const Model::real3d& v0, const Model::real3d& v1 )                                    { return Model::real3d(std::max(v0[0], v1[0]), std::max(v0[1], v1[1]), std::max(v0[2], v1[2])); }
inline Model::real3d sqr( const Model::real3d& v )                                                              { return Model::real3d(::sqr(v.c[0]), ::sqr(v.c[1]), ::sqr(v.c[2]));    }
inline Model::real   dot( const Model::real3d& v0, const Model::real3d& v1 )                                    { return v0.dot( v1 );                                                  }
inline Model::real3d cross( const Model::real3d& v0, const Model::real3d& v1 )                                  { return v0.cross( v1 );                                                }
inline Model::real3d lerp( const Model::real3d& v0, const Model::real3d& v1, Model::real64 t )                  { return (Model::real64(1.0) - t) * v0 + t * v1;                        }
inline Model::real3d lerp( const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &vt)           { return (Model::real3d(1,1,1) - vt) * v0 + vt * v1;                    }

inline Model::real3d spherical_to_cartesian( const Model::real3d& center, Model::real64 radius, Model::real64 theta, Model::real64 phi ) 
{
    Model::real64 x = std::cos(phi)*std::sin(theta);
    Model::real64 y = std::sin(phi)*std::sin(theta);
    Model::real64 z = std::cos(theta);
    return center + radius*Model::real3d(x,y,z);
}

inline Model::real4 operator +( const Model::real4& v1, const Model::real4& v2 )                                { return Model::real4(v1.c[0] + v2.c[0], v1.c[1] + v2.c[1], v1.c[2] + v2.c[2], v1.c[3] + v2.c[3] ); }
inline Model::real4 operator -( const Model::real4& v1, const Model::real4& v2 )                                { return Model::real4(v1.c[0] - v2.c[0], v1.c[1] - v2.c[1], v1.c[2] - v2.c[2], v1.c[3] - v2.c[3] ); }
inline Model::real4 operator *( const Model::real4& v1, const Model::real4& v2 )                                { return Model::real4(v1.c[0] * v2.c[0], v1.c[1] * v2.c[1], v1.c[2] * v2.c[2], v1.c[3] * v2.c[3] ); }
inline Model::real4 operator /( const Model::real4& v1, const Model::real4& v2 )                                { return Model::real4(v1.c[0] / v2.c[0], v1.c[1] / v2.c[1], v1.c[2] / v2.c[2], v1.c[3] / v2.c[3] ); }
inline Model::real4 operator *( Model::real t, const Model::real4& v )                                          { return Model::real4(t*v.c[0], t*v.c[1], t*v.c[2], t*v.c[3]);          }
inline Model::real4 operator /( Model::real4 v, Model::real t )                                                 { return Model::real4(v.c[0]/t, v.c[1]/t, v.c[2]/t, v.c[3]/t);          }
inline Model::real4 operator *( const Model::real4& v, Model::real t )                                          { return Model::real4(t*v.c[0], t*v.c[1], t*v.c[2], t*v.c[3]);          }

inline Model::real4 min( const Model::real4& v0, const Model::real4& v1 )                                       { return Model::real4(std::min(v0[0], v1[0]), std::min(v0[1], v1[1]), std::min(v0[2], v1[2]), std::min(v0[3], v1[3])); }
inline Model::real4 max( const Model::real4& v0, const Model::real4& v1 )                                       { return Model::real4(std::max(v0[0], v1[0]), std::max(v0[1], v1[1]), std::max(v0[2], v1[2]), std::max(v0[3], v1[3])); }
inline Model::real4 sqr( const Model::real4& v )                                                                { return Model::real4(::sqr(v.c[0]), ::sqr(v.c[1]), ::sqr(v.c[2]), ::sqr(v.c[3])); }
inline Model::real  dot( const Model::real4& v0, const Model::real4& v1 )                                       { return v0.dot( v1 );                                                  }
inline Model::real4 lerp( const Model::real4& v0, const Model::real4& v1, Model::real t )                       { return (Model::real(1.0) - t) * v0 + t * v1;                          }
inline Model::real4 lerp( const Model::real4 &v0, const Model::real4 &v1, const Model::real4 &vt)               { return (Model::real4(1,1,1,1) - vt) * v0 + vt * v1;                   }

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// DEBUG ASSERTS AND STRING CONVERSION UTILITIES
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
#define dprint( msg )
//#define dprint( msg ) std::cout << (msg) << "\n"
bool Model::debug = false;
#define mdout if (Model::debug) std::cout 

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 ); return false; }
#define fsc_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 ); goto error;   }
#define obj_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 ); goto error;   }
#define gph_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 ); goto error;   }
#define die_assert( bool, msg ) if ( !(bool) ) {                               std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 );               }
#define die( msg )                             {                               std::cout << "ERROR: " << msg << "\n" << std::flush; exit( 1 );               }

inline std::istream& operator >> ( std::istream& is, Model::real3& v ) 
{
    is >> v.c[0] >> v.c[1] >> v.c[2];
    return is;
}

inline std::ostream& operator << ( std::ostream& os, const Model::real4& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "," << v.c[2] << "," << v.c[3] << "]";
    return os;
}

inline std::string str( Model::_int i )
{
    return std::to_string( i );
}

inline std::string str( Model::_int64 i )
{
    return std::to_string( i );
}

inline std::string str( Model::uint u )
{
    return std::to_string( u );
}

inline std::string str( size_t s )
{
    return std::to_string( s );
}

inline std::string str( Model::real r ) 
{
    return std::to_string( r );
}

inline std::string str( Model::real64 r ) 
{
    std::ostringstream s;
    s << std::setprecision(15) << r;
    return s.str();
}

inline std::string str( Model::real3 v ) 
{
    return "[" + str(v.c[0]) + "," + str(v.c[1]) + "," + str(v.c[2]) + "]";
}

inline std::string str( Model::real3d v ) 
{
    return "[" + str(v.c[0]) + "," + str(v.c[1]) + "," + str(v.c[2]) + "]";
}

inline std::ostream& operator << ( std::ostream& os, const Model::real3& v ) 
{
    os << str( v );
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::real3d& v ) 
{
    os << str( v );
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::real2& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Matrix& m ) 
{
    os << "[ [" << m.m[0][0] << "," << m.m[0][1] << "," << m.m[0][2] << "," << m.m[0][3] << "]," << 
          "  [" << m.m[1][0] << "," << m.m[1][1] << "," << m.m[1][2] << "," << m.m[1][3] << "]," << 
          "  [" << m.m[2][0] << "," << m.m[2][1] << "," << m.m[2][2] << "," << m.m[2][3] << "]," << 
          "  [" << m.m[3][0] << "," << m.m[3][1] << "," << m.m[3][2] << "," << m.m[3][3] << "] ]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::AABB& box ) 
{
    os << box._min << ".." << box._max;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::AABBD& box ) 
{
    os << box._min << ".." << box._max;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::AABBI& box ) 
{
    os << "[" << box._min[0] << "," << box._min[1] << "," << box._min[2] << "] .. [" <<
                 box._max[0] << "," << box._max[1] << "," << box._max[2] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::AABBU64& box ) 
{
    os << "[" << box._min[0] << "," << box._min[1] << "," << box._min[2] << "] .. [" <<
                 box._max[0] << "," << box._max[1] << "," << box._max[2] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream &os, const Model::Ray &r ) 
{
    os << "{ origin=" << r.A << " direction=" << r.B << " direction_inv=" << r.B_inv << " direction_norm=" << r.B_norm << " kind=" << int(r.KIND) << 
             " solid_angle=" << r.SOLID_ANGLE << "}";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Polygon& poly ) 
{
    os << "mtl_i=" << poly.mtl_i << " vtx_cnt=" << poly.vtx_cnt << " vtx_i=" << poly.vtx_i <<
          " normal=" << poly.normal << " area=" << poly.area;
    return os;
}

inline std::string str( const Model::VolumeVoxelType type ) 
{
    switch( type ) 
    {
        case Model::VolumeVoxelType::UNKNOWN:   return "UNKNOWN";       break;
        case Model::VolumeVoxelType::FLOAT:     return "FLOAT";         break;
        case Model::VolumeVoxelType::DOUBLE:    return "DOUBLE";        break;
        case Model::VolumeVoxelType::INT16:     return "INT16";         break;
        case Model::VolumeVoxelType::INT32:     return "INT32";         break;
        case Model::VolumeVoxelType::INT64:     return "INT64";         break;
        case Model::VolumeVoxelType::VEC3F:     return "VEC3F";         break;
        case Model::VolumeVoxelType::VEC3D:     return "VEC3D";         break;
        case Model::VolumeVoxelType::MASK:      return "MASK";          break;
        case Model::VolumeVoxelType::FP16:      return "FP16";          break;
        default:                                return "<unknown>";     break;
    }
}

inline std::ostream& operator << ( std::ostream& os, const Model::VolumeVoxelType& type ) 
{
    os << str( type );
    return os;
}

inline std::string str( const Model::VolumeVoxelClass c ) 
{
    switch( c ) 
    {
        case Model::VolumeVoxelClass::UNKNOWN:          return "UNKNOWN";       break;
        case Model::VolumeVoxelClass::DENSITY:          return "DENSITY";       break;
        case Model::VolumeVoxelClass::TEMPERATURE:      return "TEMPERATURE";   break;
        case Model::VolumeVoxelClass::EMISSIVE_RGB:     return "EMISSIVE_RGB";  break;
        case Model::VolumeVoxelClass::RGB:              return "RGB";           break;
        case Model::VolumeVoxelClass::NORMAL:           return "NORMAL";        break;
        case Model::VolumeVoxelClass::IOR:              return "IOR";           break;
        case Model::VolumeVoxelClass::F0:               return "F0";            break;
        case Model::VolumeVoxelClass::G:                return "G";             break;
        case Model::VolumeVoxelClass::ATTENUATION:      return "ATTENUATION";   break;
        case Model::VolumeVoxelClass::MFPL:             return "MFPL";          break;
        default:                                        return "<unknown>";     break;
    }
}

inline std::ostream& operator << ( std::ostream& os, const Model::VolumeVoxelClass& c ) 
{
    os << str( c );
    return os;
}

inline std::string str( const Model::VolumeGridClass c ) 
{
    switch( c )
    {
        case Model::VolumeGridClass::UNKNOWN:   return "UNKNOWN";       break;
        case Model::VolumeGridClass::LEVELSET:  return "LEVELSET";      break;
        case Model::VolumeGridClass::FOGVOLUME: return "FOGVOLUME";     break;
        case Model::VolumeGridClass::STAGGERED: return "STAGGERED";     break;
        default:                                return "<unknown>";     break;
    };
}

inline std::ostream& operator << ( std::ostream& os, const Model::VolumeGridClass& c ) 
{
    os << str( c );
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::VolumeGrid& grid ) 
{
    os << "name_i=" << grid.name_i << " voxel_cnt=" << grid.voxel_cnt << " voxel_type=" << grid.voxel_type << " voxel_class=" << grid.voxel_class <<
          " grid_class=" << grid.grid_class << " world_box=" << grid.world_box << " index_box=" << grid.index_box <<
          " world_voxel_size=" << grid.world_voxel_size << " world_voxel_size_inv=" << grid.world_voxel_size_inv << 
          " world_translate=" << grid.world_translate << " voxel_i=" << grid.voxel_i;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Material& mat ) 
{
    os << "name_i="  << mat.name_i <<
          " Ka="     << mat.Ka <<
          " Kd="     << mat.Kd <<
          " Ke="     << mat.Ke <<
          " Ks="     << mat.Ks <<
          " Tf="     << mat.Tf <<
          " Tr="     << mat.Tr <<
          " Ns="     << mat.Ns <<
          " Ni="     << mat.Ni <<
          " d="      << mat.d  <<
          " illum="  << mat.illum <<
          " map_Ka_i=" << mat.map_Ka_i <<
          " map_Kd_i=" << mat.map_Kd_i <<
          " map_Ke_i=" << mat.map_Ke_i <<
          " map_Ks_i=" << mat.map_Ks_i <<
          " map_Ns_i=" << mat.map_Ns_i <<
          " map_d_i="  << mat.map_d_i <<
          " map_Bump_i=" << mat.map_Bump_i <<
          " map_refl_i=" << mat.map_refl_i;
    return os;
}

inline std::string str( const Model::GRAPH_NODE_KIND kind )
{
    return (kind == Model::GRAPH_NODE_KIND::BLOCK)  ? "BLOCK"  :
           (kind == Model::GRAPH_NODE_KIND::ASSIGN) ? "ASSIGN" : 
           (kind == Model::GRAPH_NODE_KIND::REAL)   ? "REAL"   :
           (kind == Model::GRAPH_NODE_KIND::REAL3)  ? "REAL3"  :
           (kind == Model::GRAPH_NODE_KIND::STR)    ? "STR"    :
           (kind == Model::GRAPH_NODE_KIND::ID)     ? "ID"     :
           (kind == Model::GRAPH_NODE_KIND::OP)     ? "OP"     :
           (kind == Model::GRAPH_NODE_KIND::NARY)   ? "NARY"   : "<unknown>";
}

inline std::ostream& operator << ( std::ostream& os, const Model::GRAPH_NODE_KIND& kind ) 
{
    os << str( kind );
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Texture& tex ) 
{
    os << "name_i=" << tex.name_i << " width=" << tex.width << " height=" << tex.height << " nchan=" << tex.nchan <<
                    " texel_i=" << tex.texel_i << " bump_multiplier=" << tex.bump_multiplier;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Texture::ASTC_Header& hdr ) 
{
    Model::uint w = (hdr.xsize[0] << 0) | (hdr.xsize[1] << 8) | (hdr.xsize[2] << 16);
    Model::uint h = (hdr.ysize[0] << 0) | (hdr.ysize[1] << 8) | (hdr.ysize[2] << 16);
    Model::uint d = (hdr.zsize[0] << 0) | (hdr.zsize[1] << 8) | (hdr.zsize[2] << 16);
    os << "blockdim=[" << int(hdr.blockdim_x) << "," << int(hdr.blockdim_y) << "," << int(hdr.blockdim_z) << "] size=[" << 
                          w << "," << h << "," << d << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::BVH_NODE_KIND& kind ) 
{
    os << ((kind == Model::BVH_NODE_KIND::BVH_NODE)  ? "BVH_NODE" :
           (kind == Model::BVH_NODE_KIND::POLYGON)   ? "POLYGON"  : 
           (kind == Model::BVH_NODE_KIND::INSTANCE)  ? "INSTANCE" : "<unknown>");
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::BVH_Node& bvh ) 
{
    os << "box=" << bvh.box << " left_kind=" << bvh.left_kind << " right_kind=" << bvh.right_kind << " left_i=" << bvh.left_i << " right_i=" << bvh.right_i;
    return os;
}

inline std::string Model::Graph_Node::str( const Model * model, std::string indent ) const
{
    std::ostringstream s;
    uint node_i = this - model->graph_nodes;
    s << indent << kind << " node_i=" << node_i << " parent_i=" << parent_i << " sibling_i=" << sibling_i;
    switch( kind ) 
    {
        case GRAPH_NODE_KIND::REAL:
            s << " " << u.r << "\n";
            break;

        case GRAPH_NODE_KIND::REAL3:
            s << " " << u.r3 << "\n";
            break;

        case GRAPH_NODE_KIND::STR:
            if ( u.s_i == uint(-1) ) {
                s << " <no value>\n";
            } else {
                s << " \"" << &model->strings[u.s_i] << "\"\n";
            }
            break;

        case GRAPH_NODE_KIND::ID:
        case GRAPH_NODE_KIND::OP:
            if ( u.s_i == uint(-1) ) {
                s << "<no value>\n";
            } else {
                s << " " << &model->strings[u.s_i] << "\n";
            }
            break;

        case GRAPH_NODE_KIND::BLOCK:
        case GRAPH_NODE_KIND::ASSIGN:
        case GRAPH_NODE_KIND::NARY:
        {
            s << ":\n";
            indent += "    ";
            uint this_i = this - model->graph_nodes;
            for( uint child_i = u.child_first_i; child_i != uint(-1); child_i = model->graph_nodes[child_i].sibling_i ) 
            {
                die_assert( model->graph_nodes[child_i].parent_i == this_i, "bad parent_i=" + std::to_string(model->graph_nodes[child_i].parent_i) );
                s << model->graph_nodes[child_i].str( model, indent );    
            }
            break;
        }

        default:
            break;
    }
    return s.str();
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// SYSTEM UTILITIES
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
void Model::dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ) 
{
    // in case they are not found:
    dir_name = "";
    base_name = "";
    ext_name = "";

    // ext_name 
    const int len = path.length();
    int pos = len - 1;
    for( ; pos >= 0; pos-- )
    {
        if ( path[pos] == '.' ) {
            ext_name = path.substr( pos );
            pos--;
            break;
        }
    } 
    if ( pos < 0 ) pos = len - 1;  // no ext_name, so reset for base_name
    
    // base_name
    int base_len = 0;
    for( ; pos >= 0; pos--, base_len++ )
    {
        if ( path[pos] == '/' ) {
            if ( base_len != 0 ) base_name = path.substr( pos+1, base_len );
            pos--;
            break;
        }
    }
    if ( pos < 0 ) pos = len - 1;  // no base_name, so reset for dir_name

    // dir_name is whatever's left
    if ( pos >= 0 ) dir_name = path.substr( 0, pos+1 );
}

bool Model::file_exists( std::string file_path )
{
    struct stat ss;
    return stat( file_path.c_str(), &ss ) == 0;
}

bool Model::file_read( std::string file_path, char *& start, char *& end, bool is_read_only )
{
    const char * fname = file_path.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) std::cout << "file_read() error reading " << file_path << ": " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open file " + file_path + " - open() error: " + strerror( errno ) );

    struct stat file_stat;
    int status = fstat( fd, &file_stat );
    if ( status < 0 ) {
        close( fd );
        rtn_assert( 0, "could not stat file " + std::string(fname) + " - stat() error: " + strerror( errno ) );
    }
    size_t size = file_stat.st_size;

    if ( size != 0 ) {
        // let mmap() choose an addr and make the region read-only or read/write
        int prot = PROT_READ | (is_read_only ? 0 : PROT_WRITE);
        int flags = MAP_FILE | (is_read_only ? MAP_SHARED : MAP_PRIVATE);
        void * addr = mmap( 0, size, prot, flags, fd, 0 );
        rtn_assert( addr != MAP_FAILED, "file_read() mmap() call failed - error: " + std::string( strerror( errno ) ) );
        start = reinterpret_cast<char *>( addr );
        end = start + size;
    } else {
        start = 0;
        end = start;
    }
    return true;
}

bool Model::file_read( std::string file_path, const char *& start, const char *& end )
{
    char * _start;
    char * _end;
    bool ret = file_read( file_path, _start, _end, true );
    start = _start;
    end = _end;
    return ret;
}

void Model::file_write( std::string file_path, const unsigned char * data, uint64_t byte_cnt )
{
    cmd( "rm -f '" + file_path + "'" );

    int fd = open( file_path.c_str(), O_WRONLY|O_CREAT );
    die_assert( fd >= 0, "file_write() error opening " + file_path + " for writing: " + strerror( errno ) );

    while( byte_cnt != 0 ) 
    {
        size_t _this_byte_cnt = 1024*1024*1024;
        if ( byte_cnt < _this_byte_cnt ) _this_byte_cnt = byte_cnt;
        if ( ::write( fd, data, _this_byte_cnt ) <= 0 ) {
            close( fd );
            die_assert( 0, "could not write() file " + file_path + ": " + strerror( errno ) );
        }
        byte_cnt -= _this_byte_cnt;
        data     += _this_byte_cnt;
    }
    close( fd );
    cmd( "chmod +rw " + file_path );
}

void Model::cmd( std::string c, std::string error, bool echo )
{
    if ( echo ) std::cout << c << "\n";
    if ( std::system( c.c_str() ) != 0 ) die_assert( false, "ERROR: " + error + ": " + c );
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// TOP-LEVEL METHODS FOR PARSING FILES
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
Model::Model( std::string                top_file,      
              Model::MIPMAP_FILTER       mipmap_filter, 
              Model::TEXTURE_COMPRESSION texture_compression,
              Model::BVH_TREE bvh_tree,  bool resolve_models, bool deduplicate )
{
    is_good = false;
    error_msg = "<unknown error>";
    mapped_region = nullptr;
    mapped_region_len = 0;
    fsc = nullptr;
    obj = nullptr;
    line_num = 1;

    hdr = aligned_alloc<Header>( 1 );
    memset( hdr, 0, sizeof( Header ) );
    hdr->version = VERSION;
    hdr->pos_cnt = 1;
    hdr->norm_cnt = 1;
    hdr->texcoord_cnt = 1;
    hdr->mipmap_filter = mipmap_filter;
    hdr->texture_compression = texture_compression;
    hdr->graph_node_cnt = 1;
    hdr->graph_root_i = uint(-1);
    hdr->volume_cnt = 0;
    hdr->volume_grid_cnt = 0;
    hdr->lighting_scale = 1.0;
    hdr->ambient_intensity = real3( 0.1, 0.1, 0.1 );
    hdr->sky_box_tex_i = uint(-1);
    hdr->sky_radius = 0.0;
    hdr->env_map_intensity_scale = 1.0;
    hdr->opacity_scale = 1.0;
    hdr->initial_camera_i = uint(-1);
    hdr->animation_speed = 1.0;
    hdr->background_color = real3( 0, 0, 0 );
    hdr->bvh_root_i = uint(-1);

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->obj_cnt           = 128*1024;
    max->poly_cnt          = 1024*1024;
    max->emissive_poly_cnt = 16;
    max->mtl_cnt           = 128;
    max->vtx_cnt     	   = 4*max->poly_cnt;
    max->pos_cnt     	   = max->vtx_cnt;
    max->norm_cnt    	   = max->vtx_cnt;
    max->texcoord_cnt	   = max->vtx_cnt;
    max->mipmap_filter	   = mipmap_filter;
    max->graph_node_cnt    = 1;
    max->volume_cnt        = 1;
    max->volume_grid_cnt   = 1;
    max->voxel_cnt         = 1;
    max->tex_cnt     	   = max->mtl_cnt;
    max->texel_cnt   	   = max->mtl_cnt * 128*1024;
    max->char_cnt    	   = max->obj_cnt * 128;
    max->bvh_node_cnt	   = max->poly_cnt / 2;
    max->matrix_cnt  	   = 1;
    max->inst_cnt    	   = 1;
    max->light_cnt   	   = 1;
    max->camera_cnt  	   = 1;
    max->frame_cnt   	   = 1;
    max->animation_cnt 	   = 1;

    //------------------------------------------------------------
    // Allocate arrays
    //------------------------------------------------------------
    strings           = aligned_alloc<char>(     max->char_cnt );
    objects           = aligned_alloc<Object>(   max->obj_cnt );
    polygons          = aligned_alloc<Polygon>(  max->poly_cnt );
    emissive_polygons = aligned_alloc<uint>(     max->emissive_poly_cnt );
    vertexes          = aligned_alloc<Vertex>(   max->vtx_cnt );
    positions         = aligned_alloc<real3>(    max->pos_cnt );
    normals           = aligned_alloc<real3>(    max->norm_cnt );
    texcoords         = aligned_alloc<real2>(    max->texcoord_cnt );
    materials         = aligned_alloc<Material>( max->mtl_cnt );
    graph_nodes       = aligned_alloc<Graph_Node>( max->graph_node_cnt );
    volumes           = aligned_alloc<Volume>(   max->volume_cnt );
    volume_grids      = aligned_alloc<VolumeGrid>( max->volume_grid_cnt );
    voxels            = aligned_alloc<unsigned char>( max->voxel_cnt );
    textures          = aligned_alloc<Texture>(  max->tex_cnt );
    texels            = aligned_alloc<unsigned char>( max->texel_cnt );
    bvh_nodes         = aligned_alloc<BVH_Node>( max->bvh_node_cnt );
    matrixes          = aligned_alloc<Matrix>(   max->matrix_cnt );
    instances         = aligned_alloc<Instance>( max->inst_cnt );
    lights            = aligned_alloc<Light>(    max->light_cnt );
    cameras           = aligned_alloc<Camera>(   max->camera_cnt );
    frames            = aligned_alloc<Frame>(    max->frame_cnt );
    animations        = aligned_alloc<Animation>(max->animation_cnt );

    if ( top_file == "" ) {
        //------------------------------------------------------------
        // New model.
        //------------------------------------------------------------
        is_good = true;
        return;
    }

    //------------------------------------------------------------
    // Load .fscene or .obj depending on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    std::string ext_name;
    dissect_path( top_file, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".obj" ) || ext_name == std::string( ".OBJ" ) ) {
        if ( !load_obj( top_file, dir_name ) ) return;
    } else if ( ext_name == std::string( ".fscene" ) || ext_name == std::string( ".FSCENE" ) ) {
        if ( !load_fsc( top_file, dir_name ) ) return;
    } else if ( ext_name == std::string( ".gph" ) ) {
        if ( !load_gph( top_file, dir_name, base_name, ext_name ) ) return;
    } else if ( ext_name == std::string( ".nvdb" ) ) {
        if ( !load_nvdb( top_file, dir_name, base_name, ext_name ) ) return;
    } else {
        error_msg = "unknown top file ext_name: " + ext_name;
        return;
    }

    if ( deduplicate ) {
        //------------------------------------------------------------
        // Go through each object and see if it is a duplicate of a 
        // previous object yet translated by some fixed amount.
        //------------------------------------------------------------
    }

    if ( resolve_models ) {
        //------------------------------------------------------------
        // Resolve Model Pointers in MODEL Instances
        // Do this only for unique models.
        //------------------------------------------------------------
        uint     unique_cnt = 0;
        uint *   unique_name_i = new uint[ hdr->inst_cnt ];
        Model ** uniques       = new Model*[ hdr->inst_cnt ];

        for( uint i = 0; i < hdr->inst_cnt; i++ )
        {
            Instance * instance = &instances[i];
            if ( instance->kind == INSTANCE_KIND::MODEL ) {
                uint u;
                for( u = 0; u < unique_cnt; u++ )
                {
                    if ( unique_name_i[u] == instance->model_file_name_i ) break;
                }

                Model * model;
                if ( u == unique_cnt ) {
                    // read in model for the first time
                    const char * model_file_name = &strings[instance->model_file_name_i];
                    std::string file_name = model_file_name;
                    model = new Model( file_name, false );
                    if ( !model->is_good ) {
                        error_msg = model->error_msg;
                        return;
                    }
                    unique_name_i[u] = instance->model_file_name_i;
                    uniques[u]       = model;
                    unique_cnt++;

                    // nullify any ignored_polys[]
                    auto it = model_ignored_polys.find( file_name );
                    if ( it != model_ignored_polys.end() ) {
                        std::vector<uint>& polys = *it->second;
                        for( auto pit = polys.begin(); pit != polys.end(); pit++ )
                        {
                            uint poly_i = *pit;
                            if ( poly_i < model->hdr->poly_cnt ) {
                                mdout << "Ignoring poly_i=" << poly_i << " in " << file_name << "\n";
                                model->polygons[poly_i].vtx_cnt = 0;  // makes it invisible to hit()
                            } else {
                                error_msg = "ignored poly is out of range in " + file_name;
                                return;
                            }
                        }
                    }
                } else {
                    // model already read in
                    model = uniques[u];
                }

                instance->u.model_ptr = model;
                instance->kind = INSTANCE_KIND::MODEL_PTR;

                // also transform target model's bbox to global world space and save
                // here in the instance->box.
                Matrix * M = &matrixes[instance->matrix_i];
                AABB *   t_box = &model->bvh_nodes[model->hdr->bvh_root_i].box;
                M->transform( t_box->_min, instance->box._min );
                M->transform( t_box->_max, instance->box._max );
            }
        }
    }

    //------------------------------------------------------------
    // Optionally build BVH around Instances and/or Polygons
    //------------------------------------------------------------
    if ( bvh_tree != BVH_TREE::NONE ) bvh_build( bvh_tree );

    //------------------------------------------------------------
    // Add up byte count (excluding any padding).
    //------------------------------------------------------------
    hdr->byte_cnt = uint64( 1                         ) * sizeof( hdr ) +
                    uint64( hdr->obj_cnt              ) * sizeof( objects[0] ) +
                    uint64( hdr->poly_cnt             ) * sizeof( polygons[0] ) +
                    uint64( hdr->emissive_poly_cnt    ) * sizeof( emissive_polygons[0] ) +
                    uint64( hdr->vtx_cnt              ) * sizeof( vertexes[0] ) +
                    uint64( hdr->pos_cnt              ) * sizeof( positions[0] ) +
                    uint64( hdr->norm_cnt             ) * sizeof( normals[0] ) +
                    uint64( hdr->texcoord_cnt         ) * sizeof( texcoords[0] ) +
                    uint64( hdr->mtl_cnt              ) * sizeof( materials[0] ) +
                    uint64( hdr->graph_node_cnt       ) * sizeof( graph_nodes[0] ) +
                    uint64( hdr->volume_cnt           ) * sizeof( volumes[0] ) +
                    uint64( hdr->volume_grid_cnt      ) * sizeof( volume_grids[0] ) +
                    uint64( hdr->voxel_cnt            ) * sizeof( voxels[0] ) +
                    uint64( hdr->tex_cnt              ) * sizeof( textures[0] ) +
                    uint64( hdr->texel_cnt            ) * sizeof( texels[0] ) +
                    uint64( hdr->char_cnt             ) * sizeof( strings[0] ) + 
                    uint64( hdr->bvh_node_cnt         ) * sizeof( bvh_nodes[0] ) +
                    uint64( hdr->matrix_cnt           ) * sizeof( matrixes[0] ) +
                    uint64( hdr->inst_cnt             ) * sizeof( instances[0] ) +
                    uint64( hdr->light_cnt            ) * sizeof( lights[0] ) +
                    uint64( hdr->camera_cnt           ) * sizeof( cameras[0] ) +
                    uint64( hdr->frame_cnt            ) * sizeof( frames[0] ) +
                    uint64( hdr->animation_cnt        ) * sizeof( animations[0] );

    is_good = true;
}

Model::Model( std::string model_path, bool is_compressed )
{
    is_good = false;
    mapped_region = nullptr;
    if ( !is_compressed ) {
        read_uncompressed( model_path );

    } else {
        gzFile fd = gzopen( model_path.c_str(), "r" );
        if ( fd == Z_NULL ) {
            "Could not gzopen() file " + model_path + " for reading - gzopen() error: " + strerror( errno );
            return;
        }

        //------------------------------------------------------------
        // Reader in header then individual arrays.
        //------------------------------------------------------------
        #define _read( array, type, cnt ) \
            if ( cnt == 0 ) { \
                array = nullptr; \
            } else { \
                array = aligned_alloc<type>( cnt ); \
                if ( array == nullptr ) { \
                    gzclose( fd ); \
                    error_msg = "could not allocate " #array " array"; \
                    return; \
                } \
                char * _addr = reinterpret_cast<char *>( array ); \
                for( uint64 _byte_cnt = (cnt)*sizeof(type); _byte_cnt != 0;  ) \
                { \
                    uint64 _this_byte_cnt = 1024*1024*1024; \
                    if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
                    if ( gzread( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                        gzclose( fd ); \
                        error_msg = "could not gzread() file " + model_path + " - gzread() error: " + strerror( errno ); \
                        return; \
                    } \
                    _byte_cnt -= _this_byte_cnt; \
                    _addr     += _this_byte_cnt; \
                } \
            } \

        _read( hdr,         Header,   1 );
        if ( hdr->version != VERSION ) {
            gzclose( fd );
            error_msg = "hdr->version does not match current VERSION " + std::to_string(VERSION) + ", got " + std::to_string(hdr->version) + 
                        "; please re-generate (or re-download) your .model files"; 
            die_assert( false, "aborting..." ); 
            return;
        }
        max = aligned_alloc<Header>( 1 );
        memcpy( max, hdr, sizeof( Header ) );
        _read( strings,             char,          hdr->char_cnt );
        _read( objects,             Object,        hdr->obj_cnt );
        _read( polygons,            Polygon,       hdr->poly_cnt );
        _read( emissive_polygons,   uint,          hdr->emissive_poly_cnt );
        _read( vertexes,            Vertex,        hdr->vtx_cnt );
        _read( positions,           real3,         hdr->pos_cnt );
        _read( normals,             real3,         hdr->norm_cnt );
        _read( texcoords,           real2,         hdr->texcoord_cnt );
        _read( materials,           Material,      hdr->mtl_cnt );
        _read( graph_nodes,         Graph_Node,    hdr->graph_node_cnt );
        _read( volumes,             Volume,        hdr->volume_cnt );
        _read( volume_grids,        VolumeGrid,    hdr->volume_grid_cnt );
        _read( voxels,              unsigned char, hdr->voxel_cnt );
        _read( textures,            Texture,       hdr->tex_cnt );
        _read( texels,              unsigned char, hdr->texel_cnt );
        _read( bvh_nodes,           BVH_Node,      hdr->bvh_node_cnt );
        _read( matrixes,            Matrix,        hdr->matrix_cnt );
        _read( instances,           Instance,      hdr->inst_cnt );
        _read( lights,              Light,         hdr->light_cnt );
        _read( cameras,             Camera,        hdr->camera_cnt );
        _read( frames,              Frame,         hdr->frame_cnt );
        _read( animations,          Animation,     hdr->animation_cnt );

        gzclose( fd );
    }

    //------------------------------------------------------------
    // Recreate maps.
    //------------------------------------------------------------
    for( uint i = 0; i < hdr->obj_cnt; i++ )
    {
        std::string name = &strings[objects[i].name_i];
        name_to_obj_i[name] = i;
        mdout << "obj_i=" << i << " name_i=" << objects[i].name_i << " name=" << name << "\n";
    }

    for( uint i = 0; i < hdr->mtl_cnt; i++ )
    {
        std::string name = &strings[materials[i].name_i];
        name_to_mtl_i[name] = i;
        mdout << "mtl_i=" << i << " name_i=" << materials[i].name_i << " name=" << name << "\n";
    }

    for( uint i = 0; i < hdr->tex_cnt; i++ )
    {
        std::string name = &strings[textures[i].name_i];
        name_to_tex_i[name] = i;
        mdout << "tex_i=" << i << " name_i=" << textures[i].name_i << " name=" << name << "\n";
    }

    volume_to_emissive = (hdr->volume_cnt == 0) ? nullptr : (new const uint64_t *[hdr->volume_cnt]);
    for( uint i = 0; i < hdr->volume_cnt; i++ )
    {
        volume_to_emissive[i] = nullptr;

        //------------------------------------------------------------
        // Precompute emissive information if this is an emissive volume.
        // We store bitmasks for 4x4x4 supervoxels in uint64_t.
        //------------------------------------------------------------
        const Volume * volume = &volumes[i];
        uint t_grid_i = volume->voxel_class_grid_i( this, Model::VolumeVoxelClass::TEMPERATURE );
        if ( t_grid_i == uint(-1) ) continue;
        uint d_grid_i = volume->voxel_class_grid_i( this, Model::VolumeVoxelClass::DENSITY );
        if ( d_grid_i == uint(-1) ) continue;
        const VolumeGrid * d_grid_ptr = &volume_grids[d_grid_i];
        const VolumeGrid * t_grid_ptr = &volume_grids[t_grid_i];

        _int x_min = d_grid_ptr->index_box._min[0];
        _int y_min = d_grid_ptr->index_box._min[0];
        _int z_min = d_grid_ptr->index_box._min[0];
        _int x_max = d_grid_ptr->index_box._max[0];
        _int y_max = d_grid_ptr->index_box._max[0];
        _int z_max = d_grid_ptr->index_box._max[0];
        uint x_cnt = x_max - x_min + 1;
        uint y_cnt = y_max - y_min + 1;
        uint z_cnt = z_max - z_min + 1;
        uint super_x_cnt = (x_cnt + 3) / 4;
        uint super_y_cnt = (y_cnt + 3) / 4;
        uint super_xy_cnt = super_x_cnt * super_y_cnt;
        uint super_z_cnt = (z_cnt + 3) / 4;
        uint super_cnt = super_x_cnt * super_y_cnt * super_z_cnt;
        uint64_t * super_voxels = new uint64_t[super_cnt];
        memset( super_voxels, 0, super_cnt * sizeof(super_voxels[0]) );
        volume_to_emissive[i] = super_voxels;

        //------------------------------------------------------------
        // Now set bits for each voxel with density > 0 and temperature > 0.
        //------------------------------------------------------------
        for( _int z = z_min; z <= z_max; z++ )
        {
            for( _int y = y_min; y <= y_max; y++ )
            {
                for( _int x = x_min; x <= x_max; x++ )
                {
                    if ( t_grid_ptr->real_value( this, x, y, z ) > 0.0 && d_grid_ptr->real_value( this, x, y, z ) > 0.0 ) {
                        uint super_x = (x - x_min) / 4;
                        uint super_y = (y - y_min) / 4;
                        uint super_z = (z - z_min) / 4;
                        uint x_off   = (x - x_min) % 4;
                        uint y_off   = (y - y_min) % 4;
                        uint z_off   = (z - x_min) % 4;
                        uint super_i = super_z*super_xy_cnt + super_y*super_x_cnt + super_x;
                        uint bit_i   = z_off*16 + y_off*4 + x_off;
                        super_voxels[super_i] |= 1ULL << bit_i;
                    }
                }
            }
        }
    }

    is_good = true;
}

Model::~Model()
{
    if ( mapped_region != nullptr ) {
        munmap( mapped_region, mapped_region_len );
        mapped_region = nullptr;
    } else {
        delete strings;
        delete objects;
        delete polygons;
        delete emissive_polygons;
        delete vertexes;
        delete positions;
        delete normals;
        delete texcoords;
        delete materials;
        delete graph_nodes;
        delete volumes;
        delete volume_grids;
        delete voxels;
        delete textures;
        delete texels;
        delete bvh_nodes;
        delete matrixes;
        delete instances;
        delete lights;
        delete cameras;
        delete frames;
        delete animations;
    }
}

bool Model::write( std::string model_path, bool is_compressed ) 
{
    if ( !is_compressed ) return write_uncompressed( model_path );

    gzFile fd = gzopen( model_path.c_str(), "w" );
    rtn_assert( fd != Z_NULL, "could not gzopen() file " + model_path + " for writing - gzopen() error: " + strerror( errno ) );

    //------------------------------------------------------------
    // Write out header than individual arrays.
    //------------------------------------------------------------
    #define _write( addr, byte_cnt ) \
    { \
        char * _addr = reinterpret_cast<char *>( addr ); \
        for( uint64 _byte_cnt = byte_cnt; _byte_cnt != 0;  ) \
        { \
            uint64 _this_byte_cnt = 1024*1024*1024; \
            if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
            if ( gzwrite( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                gzclose( fd ); \
                rtn_assert( 0, "could not gzwrite() file " + model_path + " - gzwrite() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _write( hdr,                1                       * sizeof(hdr[0]) );
    _write( strings,            hdr->char_cnt           * sizeof(strings[0]) );
    _write( objects,            hdr->obj_cnt            * sizeof(objects[0]) );
    _write( polygons,           hdr->poly_cnt           * sizeof(polygons[0]) );
    _write( emissive_polygons,  hdr->emissive_poly_cnt  * sizeof(emissive_polygons[0]) );
    _write( vertexes,           hdr->vtx_cnt            * sizeof(vertexes[0]) );
    _write( positions,          hdr->pos_cnt            * sizeof(positions[0]) );
    _write( normals,            hdr->norm_cnt           * sizeof(normals[0]) );
    _write( texcoords,          hdr->texcoord_cnt       * sizeof(texcoords[0]) );
    _write( materials,   	hdr->mtl_cnt       	* sizeof(materials[0]) );
    _write( graph_nodes,   	hdr->graph_node_cnt  	* sizeof(graph_nodes[0]) );
    _write( volumes,            hdr->volume_cnt         * sizeof(volumes[0]) );
    _write( volume_grids,       hdr->volume_grid_cnt    * sizeof(volume_grids[0]) );
    _write( voxels,             hdr->voxel_cnt          * sizeof(voxels[0]) );
    _write( textures,    	hdr->tex_cnt       	* sizeof(textures[0]) );
    _write( texels,      	hdr->texel_cnt     	* sizeof(texels[0]) );
    _write( bvh_nodes,   	hdr->bvh_node_cnt  	* sizeof(bvh_nodes[0]) );
    _write( matrixes,    	hdr->matrix_cnt    	* sizeof(matrixes[0]) );
    _write( instances,   	hdr->inst_cnt      	* sizeof(instances[0]) );
    _write( lights,      	hdr->light_cnt     	* sizeof(lights[0]) );
    _write( cameras,     	hdr->camera_cnt    	* sizeof(cameras[0]) );
    _write( frames,      	hdr->frame_cnt     	* sizeof(frames[0]) );
    _write( animations,  	hdr->animation_cnt 	* sizeof(animations[0]) );

    gzclose( fd );
    return true;
}

bool Model::replace_materials( std::string mtl_file_path )
{
    return load_mtl( mtl_file_path, "", true );
}

// returns array of T on a page boundary
template<typename T>
T * Model::aligned_alloc( Model::uint64 cnt )
{
    void * mem = nullptr;
    posix_memalign( &mem, getpagesize(), cnt*sizeof(T) );
    return reinterpret_cast<T *>( mem );
}

// reallocate array if we are about to exceed its current size
template<typename T>
inline void Model::perhaps_realloc( T *& array, const Model::uint64& hdr_cnt, Model::uint64& max_cnt, Model::uint64 add_cnt )
{
    while( (hdr_cnt + add_cnt) > max_cnt ) {
        void * mem = nullptr;
        uint64 old_max_cnt = max_cnt;
        max_cnt *= 2;
        if ( max_cnt < old_max_cnt ) {
            max_cnt = uint(-1);
        }
        posix_memalign( &mem, getpagesize(), max_cnt*sizeof(T) );
        memcpy( mem, array, hdr_cnt*sizeof(T) );
        delete array;
        array = reinterpret_cast<T *>( mem );
    }
}

inline uint Model::make_string( std::string s )
{
    uint s_len = s.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    uint s_i = hdr->char_cnt;
    char * to_s = &strings[hdr->char_cnt];
    hdr->char_cnt += s_len + 1;
    memcpy( to_s, s.c_str(), s_len+1 );
    return s_i;
}

inline uint Model::make_texture( std::string name, std::string dir )
{
    Texture * texture;
    if ( !load_tex( name.c_str(), dir, texture ) ) return uint(-1);
    return texture - textures;
}

inline uint Model::make_texture( std::string name, const unsigned char * data, uint nx, uint ny, uint nchan )
{
    Texture * texture;
    if ( !load_tex( name.c_str(), "", texture, data, nx, ny, nchan ) ) return uint(-1);
    return texture - textures;
}

inline uint Model::make_constant_texture( std::string name, const real4& value )
{
    unsigned char data[4];
    for( uint c = 0; c < 4; c++ )
    {
        real rgba = value[c];
        data[c] = (rgba >= 1.0) ? 255 : (rgba < 0.0) ? 0 : uint(rgba*255);
    }
    return make_texture( name.c_str(), data, 1, 1, 4 );
}

inline uint Model::make_constant_texture( std::string name, const real3& value )
{
    return make_constant_texture( name, Model::real4(value[0], value[1], value[2], 1.0) ); 
}

inline uint Model::make_material( std::string name, real3 Ka, real3 Kd, real3 Ke, real3 Ks, real3 Tf, real Tr, real Ns, real Ni, real d, real illum,
                                  uint map_Ka_i, uint map_Kd_i, uint map_Ke_i, uint map_Ks_i, uint map_Ns_i, uint map_d_i, uint map_Bump_i, uint map_refl_i )
{
    perhaps_realloc<Material>( materials, hdr->mtl_cnt, max->mtl_cnt, 1 );
    uint mtl_i = hdr->mtl_cnt++;
    Material * material = &materials[mtl_i];
    material->name_i = make_string( name );
    material->Ka = Ka;
    material->Kd = Kd;
    material->Ke = Ke;
    material->Ks = Ks;
    material->Tf = Tf;
    material->Tr = Tr;
    material->Ns = Ns;
    material->Ni = Ni;
    material->d = d;
    material->illum = illum;
    material->map_Ka_i = map_Ka_i;
    material->map_Kd_i = map_Kd_i;
    material->map_Ke_i = map_Ke_i;
    material->map_Ks_i = map_Ks_i;
    material->map_Ns_i = map_Ns_i;
    material->map_d_i = map_d_i;
    material->map_Bump_i = map_Bump_i;
    material->map_refl_i = map_refl_i;
    return mtl_i;
}

inline uint Model::make_emissive_material( std::string name, real3 Ke, uint map_Ke_i, real d )
{
    real3 zero3( 0.0, 0.0, 0.0 );
    uint  no_map = uint(-1);
    return make_material( name, zero3, zero3, Ke, zero3, zero3, 0.0, 0.0, 0.0, d, 0.0, no_map, no_map, map_Ke_i, no_map, no_map, no_map, no_map, no_map );                         
}

inline uint Model::make_diffuse_material( std::string name, real3 Kd, uint map_Kd_i, real d, uint map_d_i ) 
{
    real3 zero3( 0.0, 0.0, 0.0 );
    uint  no_map = uint(-1);
    return make_material( name, zero3, Kd, zero3, zero3, zero3, 0.0, 0.0, 0.0, d, 0.0, no_map, map_Kd_i, no_map, no_map, no_map, map_d_i, no_map, no_map );                         
}

inline uint Model::make_vertex( const real3& p, const real3& n, const real2& uv ) 
{
    perhaps_realloc<real3>( positions, hdr->pos_cnt, max->pos_cnt, 1 );
    perhaps_realloc<real3>( normals, hdr->norm_cnt, max->norm_cnt, 1 );
    perhaps_realloc<real2>( texcoords, hdr->texcoord_cnt, max->texcoord_cnt, 1 );
    perhaps_realloc<Vertex>( vertexes, hdr->vtx_cnt, max->vtx_cnt, 1 );
    uint vtx_i = hdr->vtx_cnt++;
    Vertex * vertex = &vertexes[vtx_i];
    vertex->v_i = hdr->pos_cnt;
    vertex->vn_i = hdr->norm_cnt;
    vertex->vt_i = hdr->texcoord_cnt;
    positions[hdr->pos_cnt++] = p;
    normals[hdr->norm_cnt++] = n;
    texcoords[hdr->texcoord_cnt++] = uv;
    mdout << "Model::make_vertex:     vtx_i =" << vtx_i << " p=" << p << " n=" << n << " uv=" << uv << "\n";
    return vtx_i;
}

inline uint Model::make_polygon( uint vtx_cnt, const real3 p[], const real3 n[], const real2 uv[], uint mtl_i )
{
    perhaps_realloc<Polygon>( polygons, hdr->poly_cnt, max->poly_cnt, 1 );
    uint poly_i = hdr->poly_cnt++;
    Polygon * poly = &polygons[poly_i];
    poly->vtx_cnt = vtx_cnt;
    poly->vtx_i = hdr->vtx_cnt;
    poly->mtl_i = mtl_i;
    poly->normal = (p[1] - p[0]).cross( p[2] - p[0] );
    real len = poly->normal.length();
    poly->area = ::divby2(len);                       // correct only for triangles
    poly->normal /= len;
    mdout << "Model::make_polygon: poly_i=" << poly_i << " mtl_i=" << mtl_i << " normal=" << poly->normal << " area=" << poly->area << "\n";
    for( uint i = 0; i < vtx_cnt; i++ )
    {
        mdout << "    vtx" << i << " p=" << p[i] << " n=" << n[i] << " uv=" << uv[i] << "\n";
        make_vertex( p[i], n[i], uv[i] );
    }
    if ( mtl_i != uint(-1) && materials[mtl_i].is_emissive() ) {
        perhaps_realloc<uint>( emissive_polygons, hdr->emissive_poly_cnt, max->emissive_poly_cnt, 1 );
        emissive_polygons[hdr->emissive_poly_cnt++] = poly_i;
    }
    return poly_i;
}

inline uint Model::make_polygon( uint vtx_cnt, const real3 p[], uint mtl_i )
{
    perhaps_realloc<Polygon>( polygons, hdr->poly_cnt, max->poly_cnt, 1 );
    uint poly_i = hdr->poly_cnt++;
    Polygon * poly = &polygons[poly_i];
    poly->vtx_cnt = vtx_cnt;
    poly->vtx_i = hdr->vtx_cnt;
    poly->mtl_i = mtl_i;
    poly->normal = (p[1] - p[0]).cross( p[2] - p[0] );
    real len = poly->normal.length();
    poly->area = ::divby2(len);                       // correct only for triangles
    poly->normal /= len;
    mdout << "Model::make_polygon: poly_i=" << poly_i << " mtl_i=" << mtl_i << " normal=" << poly->normal << " area=" << poly->area << "\n";
    for( uint i = 0; i < vtx_cnt; i++ )
    {
        real u = (i == 1) ? 1.0 : 0.0;
        real v = (i == 2) ? 1.0 : 0.0;
        real2 uv( u, v );
        mdout << "    vtx" << i << " p=" << p[i] << " uv=" << uv << "\n";
        make_vertex( p[i], poly->normal, uv );
    }
    if ( mtl_i != uint(-1) && materials[mtl_i].is_emissive() ) {
        perhaps_realloc<uint>( emissive_polygons, hdr->emissive_poly_cnt, max->emissive_poly_cnt, 1 );
        emissive_polygons[hdr->emissive_poly_cnt++] = poly_i;
    }
    return poly_i;
}

inline uint Model::make_sphere( const real3& center, real radius, uint nlong, uint nlat, uint mtl_i ) 
{
    real delta_u = rcp(real(nlong-1));
    real delta_v = rcp(real(nlat-1));
    uint first_poly_i = hdr->poly_cnt;
    for ( uint i = 0; i < nlong; i++ ) {
        for ( uint j = 0; j < nlat; j++ ) {
            real u = real(i)/real(nlong-1);
            real v = real(j)/real(nlat-1);
            real phi = PI2*u;
            real theta = PI*v;
            real next_phi = phi + PI2*delta_u;
            real next_theta = theta + PI*delta_v;
            real3 p00 = spherical_to_cartesian( center, radius, theta, phi );
            real3 p01 = spherical_to_cartesian( center, radius, next_theta, phi );
            real3 p10 = spherical_to_cartesian( center, radius, theta, next_phi );
            real3 p11 = spherical_to_cartesian( center, radius, next_theta, next_phi );
            real3 n00 = (p00-center).normalized();
            real3 n01 = (p01-center).normalized();
            real3 n10 = (p10-center).normalized();
            real3 n11 = (p11-center).normalized();
            real3 p1[3]  = { p00, p10, p01 };
            real3 p2[3]  = { p11, p10, p01 };
            real3 n1[3]  = { n00, n10, n01 };
            real3 n2[3]  = { n11, n10, n01 };
            real2 uv00   = real2(u, v);
            real2 uv01   = real2(u, v+delta_v);
            real2 uv10   = real2(u+delta_u, v);
            real2 uv11   = real2(u+delta_u, v+delta_v);
            real2 uv1[3] = { uv00, uv10, uv01 }; 
            real2 uv2[3] = { uv11, uv10, uv01 };
            make_polygon( 3, p1, n1, uv1, mtl_i );
            make_polygon( 3, p2, n2, uv2, mtl_i );
        }
    }
    return first_poly_i;
}

inline uint Model::make_box( const real3 front[4], const real3 back[4], uint mtl_i )
{
    //------------------------------------------------------------
    // We assume that front[0] aligns with back[0] as the bottom-left for each (as viewed from the front).
    // Use make_polygon() to create triangles for each face.
    // Be careful about face vertex order to ensure faces point outward.
    //------------------------------------------------------------
    uint first_poly_i = hdr->poly_cnt;

    const real3 tri0[3]  = { front[0], front[1], front[2] }; make_polygon( 3, tri0,  mtl_i );    // front
    const real3 tri1[3]  = { front[0], front[2], front[3] }; make_polygon( 3, tri1,  mtl_i );    // front
    const real3 tri2[3]  = { back[1],  back[0],  back[3]  }; make_polygon( 3, tri2,  mtl_i );    // back  
    const real3 tri3[3]  = { back[1],  back[3],  back[2]  }; make_polygon( 3, tri3,  mtl_i );    // back  
    const real3 tri4[3]  = { front[3], front[2], back[2]  }; make_polygon( 3, tri4,  mtl_i );    // top
    const real3 tri5[3]  = { front[3], back[2],  back[3]  }; make_polygon( 3, tri5,  mtl_i );    // top
    const real3 tri6[3]  = { back[0],  back[1],  front[1] }; make_polygon( 3, tri6,  mtl_i );    // bottom
    const real3 tri7[3]  = { back[0],  front[1], front[0] }; make_polygon( 3, tri7,  mtl_i );    // bottom 
    const real3 tri8[3]  = { back[0],  front[0], front[3] }; make_polygon( 3, tri8,  mtl_i );    // left
    const real3 tri9[3]  = { back[0],  front[3], back[3]  }; make_polygon( 3, tri9,  mtl_i );    // left
    const real3 tri10[3] = { front[1], back[1],  back[2]  }; make_polygon( 3, tri10, mtl_i );    // right
    const real3 tri11[3] = { front[1], back[2],  front[2] }; make_polygon( 3, tri11, mtl_i );    // right

    return first_poly_i;
}

inline uint Model::make_aabox( const AABB& box, uint mtl_i )
{
    //------------------------------------------------------------
    // Use make_box() after we set up front and back faces
    // the way it wants them.
    //------------------------------------------------------------
    const real3& min = box._min;
    const real3& max = box._max;
    real3 front[4] = { real3( min[0], min[1], min[2] ),
                       real3( max[0], min[1], min[2] ),
                       real3( max[0], max[1], min[2] ),
                       real3( min[0], max[1], min[2] ) };
    real3 back[4]  = { real3( min[0], min[1], max[2] ),
                       real3( max[0], min[1], max[2] ),
                       real3( max[0], max[1], max[2] ),
                       real3( min[0], max[1], max[2] ) };
    return make_box( front, back, mtl_i );
}

inline uint Model::make_skybox( const real3& center, real radius, const real3& up, uint nlong, uint nlat, uint mtl_i ) 
{
    assert(center == real3(0,0,0)); // temporary requirement until this code gets merged with previous routine

    real3 b1, b0;
    branchlessONB(up, b1, b0) ;

    uint first_poly_i = hdr->poly_cnt;
    for ( uint i = 0; i < nlong; i++ ) {
        for ( uint j = 0; j < nlat; j++ ) {
            real x0, x1, x2, x3;
            real y0, y1, y2, y3;
            real z0, z1, z2, z3;
            real A1, A2, A3, A4;
            real u0, v0, u1, v1;
            u0 = real(i) / real(nlong);
            v0 = real(j) / real(nlat);
            u1 = real(i + 1) / real(nlong);
            v1 = real(j + 1) / real(nlat);
            concentric_square_disk(u0, v0, x0, y0);
            concentric_square_disk(u1, v0, x1, y1);
            concentric_square_disk(u1, v1, x2, y2);
            concentric_square_disk(u0, v1, x3, y3);

            A1 = area_2D_triangle(x0, y0, x1, y1, x2, y2);
            A2 = area_2D_triangle(x0, y0, x2, y2, x3, y3);
            A3 = area_2D_triangle(x0, y0, x1, y1, x3, y3);
            A4 = area_2D_triangle(x1, y1, x2, y2, x3, y3);
            z0 = x0*x0 + y0*y0;
            z1 = x1*x1 + y1*y1;
            z2 = x2*x2 + y2*y2;
            z3 = x3*x3 + y3*y3;
            if (z0 >= 0.99)
                z0 = 0.0;
            else
                z0 = std::sqrt(1.0-z0);
            if (z1 >= 0.99)
                z1 = 0.0;
            else
                z1 = std::sqrt(1.0-z1);
            if (z2 >= 0.99)
                z2 = 0.0;
            else
                z2 = std::sqrt(1.0-z2);
            if (z3 >= 0.99)
                z3 = 0.0;
            else
                z3 = std::sqrt(1.0-z3);
            real3 p0(x0,  y0, z0);
            real3 p1(x1,  y1, z1);
            real3 p2(x2,  y2, z2);
            real3 p3(x3,  y3, z3);
            p0 = p0[0]*b0 + p0[1]*b1 + p0[2]*up;
            p1 = p1[0]*b0 + p1[1]*b1 + p1[2]*up;
            p2 = p2[0]*b0 + p2[1]*b1 + p2[2]*up;
            p3 = p3[0]*b0 + p3[1]*b1 + p3[2]*up;
            p0 *= radius;
            p1 *= radius;
            p2 *= radius;
            p3 *= radius;
            real3 n0 = (center-p0).normalized();
            real3 n1 = (center-p1).normalized();
            real3 n2 = (center-p2).normalized();
            real3 n3 = (center-p3).normalized();
            real2 uv0 = real2(u0, v0);
            real2 uv1 = real2(u1, v0);
            real2 uv2 = real2(u1, v1);
            real2 uv3 = real2(u0, v1);
            real3 p[3];
            real3 n[3];
            real2 uv[3];
            if ( A1/A2 + A2/A1 < A3/A4 + A4/A3 ) {
                p[0]  = p0;  p[1]  = p1;  p[2]  = p2; 
                n[0]  = n0;  n[1]  = n1;  n[2]  = n2; 
                uv[0] = uv0; uv[1] = uv1; uv[2] = uv2;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = p0;  p[1]  = p2;  p[2]  = p3; 
                n[0]  = n0;  n[1]  = n2;  n[2]  = n3; 
                uv[0] = uv0; uv[1] = uv2; uv[2] = uv3;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = -p0;  p[1] = -p2; p[2]  = -p1; 
                n[0]  = -n0;  n[1] = -n2; n[2]  = -n1; 
                uv[0] = uv0; uv[1] = uv2; uv[2] = uv1;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = -p0;  p[1] = -p3; p[2]  = -p2; 
                n[0]  = -n0;  n[1] = -n3; n[2]  = -n2; 
                uv[0] = uv0; uv[1] = uv3; uv[2] = uv2;
                make_polygon( 3, p, n, uv, mtl_i );
            } else {
                p[0]  = p0;  p[1]  = p1;  p[2]  = p3; 
                n[0]  = n0;  n[1]  = n1;  n[2]  = n3; 
                uv[0] = uv0; uv[1] = uv1; uv[2] = uv3;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = p1;  p[1]  = p2;  p[2]  = p3; 
                n[0]  = n1;  n[1]  = n2;  n[2]  = n3; 
                uv[0] = uv1; uv[1] = uv2; uv[2] = uv3;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = -p0;  p[1] = -p3; p[2]  = -p1; 
                n[0]  = -n0;  n[1] = -n3; n[2]  = -n1; 
                uv[0] = uv0; uv[1] = uv3; uv[2] = uv1;
                make_polygon( 3, p, n, uv, mtl_i );

                p[0]  = -p1;  p[1] = -p3; p[2]  = -p2; 
                n[0]  = -n1;  n[1] = -n3; n[2]  = -n2; 
                uv[0] = uv1; uv[1] = uv3; uv[2] = uv2;
                make_polygon( 3, p, n, uv, mtl_i );
            }
        }
    }
    return first_poly_i;
}

inline uint Model::make_matrix( void )
{
    perhaps_realloc<Matrix>( matrixes, hdr->matrix_cnt, max->matrix_cnt, 1 );
    uint matrix_i = hdr->matrix_cnt++;
    matrixes[matrix_i].identity();
    return matrix_i;
}

inline uint Model::make_instance( std::string name, std::string model_name, Model * model, uint matrix_i, std::string file_name )
{
    perhaps_realloc<Instance>( instances, hdr->inst_cnt, max->inst_cnt, 1 );
    uint inst_i = hdr->inst_cnt++;
    Instance * instance = &instances[inst_i];
    instance->kind = INSTANCE_KIND::MODEL_PTR;
    instance->name_i = make_string( name );
    instance->model_name_i = make_string( model_name );
    instance->model_file_name_i = make_string( file_name );
    instance->u.model_ptr = model;
    if ( matrix_i == uint(-1) ) matrix_i = make_matrix();
    instance->matrix_i = matrix_i;
    instance->matrix_inv_i = make_matrix();
    instance->matrix_inv_trans_i = make_matrix();
    matrixes[instance->matrix_i].invert( matrixes[instance->matrix_inv_i] );
    matrixes[instance->matrix_inv_i].transpose( matrixes[instance->matrix_inv_trans_i] );
    if ( model->hdr->bvh_root_i != uint(-1) ) {
        AABB * t_box = &model->bvh_nodes[model->hdr->bvh_root_i].box;
        matrixes[instance->matrix_i].transform( t_box->_min, instance->box._min );
        matrixes[instance->matrix_i].transform( t_box->_max, instance->box._max );
    }
    return inst_i;
}

bool Model::write_uncompressed( std::string model_path ) 
{
    int fd = open( model_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( fd < 0 ) std::cout << "open() for write error: " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open() file " + model_path + " for writing - open() error: " + strerror( errno ) );

    //------------------------------------------------------------
    // Write out header than individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    size_t page_size = getpagesize();

    #define _uwrite( addr, byte_cnt ) \
    { \
        size_t _byte_cnt = byte_cnt; \
        _byte_cnt += _byte_cnt % page_size; \
        char * _addr = reinterpret_cast<char *>( addr ); \
        for( ; _byte_cnt != 0;  ) \
        { \
            uint _this_byte_cnt = 1024*1024*1024; \
            if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
            if ( ::write( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                close( fd ); \
                rtn_assert( 0, "could not write() file " + model_path + " - write() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _uwrite( hdr,               1                       * sizeof(hdr[0]) );
    _uwrite( strings,           hdr->char_cnt           * sizeof(strings[0]) );
    _uwrite( objects,           hdr->obj_cnt            * sizeof(objects[0]) );
    _uwrite( polygons,          hdr->poly_cnt           * sizeof(polygons[0]) );
    _uwrite( emissive_polygons, hdr->emissive_poly_cnt  * sizeof(emissive_polygons[0]) );
    _uwrite( vertexes,          hdr->vtx_cnt            * sizeof(vertexes[0]) );
    _uwrite( positions,         hdr->pos_cnt            * sizeof(positions[0]) );
    _uwrite( normals,           hdr->norm_cnt           * sizeof(normals[0]) );
    _uwrite( texcoords,         hdr->texcoord_cnt       * sizeof(texcoords[0]) );
    _uwrite( materials,   	hdr->mtl_cnt       	* sizeof(materials[0]) );
    _uwrite( graph_nodes,   	hdr->graph_node_cnt  	* sizeof(graph_nodes[0]) );
    _uwrite( volumes,           hdr->volume_cnt         * sizeof(volumes[0]) );
    _uwrite( volume_grids,      hdr->volume_grid_cnt    * sizeof(volume_grids[0]) );
    _uwrite( voxels,            hdr->voxel_cnt          * sizeof(voxels[0]) );
    _uwrite( textures,    	hdr->tex_cnt       	* sizeof(textures[0]) );
    _uwrite( texels,      	hdr->texel_cnt     	* sizeof(texels[0]) );
    _uwrite( bvh_nodes,   	hdr->bvh_node_cnt  	* sizeof(bvh_nodes[0]) );
    _uwrite( matrixes,    	hdr->matrix_cnt    	* sizeof(matrixes[0]) );
    _uwrite( instances,   	hdr->inst_cnt      	* sizeof(instances[0]) );
    _uwrite( lights,      	hdr->light_cnt     	* sizeof(lights[0]) );
    _uwrite( cameras,     	hdr->camera_cnt    	* sizeof(cameras[0]) );
    _uwrite( frames,      	hdr->frame_cnt     	* sizeof(frames[0]) );
    _uwrite( animations,  	hdr->animation_cnt 	* sizeof(animations[0]) );

    fsync( fd ); // flush
    close( fd );
    std::string c = "chmod 666 ";
    c += model_path.c_str();
    cmd( c );

    if ( 0 ) {
        for( uint32_t i = 0; i < hdr->bvh_node_cnt; i++ )
        {
            std::cout << "bvh[" << i << "]: " << bvh_nodes[i] << "\n";
        }
    }
    return true;
}

bool Model::read_uncompressed( std::string model_path )
{
    char * start;
    char * end;
    if ( !file_read( model_path, start, end, false ) ) return false; // still some code that modifies data, so need read/write
    mapped_region = start;
    mapped_region_len = end - start;

    //------------------------------------------------------------
    // Write out header than individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    char * _addr = start;
    size_t page_size = getpagesize();

    #define _uread( array, type, cnt ) \
        if ( (cnt) == 0 ) { \
            array = nullptr; \
        } else { \
            array = reinterpret_cast<type *>( _addr ); \
            size_t _byte_cnt = (cnt)*sizeof(type); \
            _byte_cnt += _byte_cnt % page_size; \
            _addr += _byte_cnt; \
        } \

    _uread( hdr,         Header,   1 );
    if ( hdr->version != VERSION ) {
        rtn_assert( 0, "hdr->version does not match current VERSION " + std::to_string(VERSION) + ", got " + 
                       std::to_string(hdr->version) + "; please re-generate (or re-download) your .model files" );  
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _uread( strings,             char,          hdr->char_cnt );
    _uread( objects,             Object,        hdr->obj_cnt );
    _uread( polygons,            Polygon,       hdr->poly_cnt );
    _uread( emissive_polygons,   uint,          hdr->emissive_poly_cnt );
    _uread( vertexes,            Vertex,        hdr->vtx_cnt );
    _uread( positions,           real3,         hdr->pos_cnt );
    _uread( normals,             real3,         hdr->norm_cnt );
    _uread( texcoords,           real2,         hdr->texcoord_cnt );
    _uread( materials,           Material,      hdr->mtl_cnt );
    _uread( graph_nodes,         Graph_Node,    hdr->graph_node_cnt );
    _uread( volumes,             Volume,        hdr->volume_cnt );
    _uread( volume_grids,        VolumeGrid,    hdr->volume_grid_cnt );
    _uread( voxels,              unsigned char, hdr->voxel_cnt );
    _uread( textures,            Texture,       hdr->tex_cnt );
    _uread( texels,              unsigned char, hdr->texel_cnt );
    _uread( bvh_nodes,           BVH_Node,      hdr->bvh_node_cnt );
    _uread( matrixes,            Matrix,        hdr->matrix_cnt );
    _uread( instances,           Instance,      hdr->inst_cnt );
    _uread( lights,              Light,         hdr->light_cnt );
    _uread( cameras,             Camera,        hdr->camera_cnt );
    _uread( frames,              Frame,         hdr->frame_cnt );
    _uread( animations,          Animation,     hdr->animation_cnt );

    is_good = true;

    return true;
}

bool Model::load_fsc( std::string fsc_file, std::string dir_name )
{
    (void)dir_name;

    //------------------------------------------------------------
    // Map in .fscene file
    //------------------------------------------------------------
    line_num = 1;
    fsc = nullptr;
    uint initial_camera_name_i = uint(-1);
    if ( !file_read( fsc_file, fsc_start, fsc_end ) ) goto error;
    fsc = fsc_start;

    //------------------------------------------------------------
    // Parse top dictionary.
    //------------------------------------------------------------
    hdr->tone_key = 0.2;
    hdr->tone_white = 3.0;
    hdr->tone_avg_luminance = 0.0;
    hdr->force_tone_none = false;
    hdr->force_tone_avg_luminance = false;
    hdr->tex_specular_is_orm = false;
    if ( !expect_char( '{', fsc, fsc_end, true ) ) goto error;
    for( ;; )
    {
        skip_whitespace( fsc, fsc_end );
        fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
        if ( *fsc == '}' ) {
            if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
            skip_whitespace( fsc, fsc_end );
            fsc_assert( fsc == fsc_end, "unexpected characters at end of " + fsc_file );
            if ( initial_camera_name_i == uint(-1) ) {
                hdr->initial_camera_i = uint(-1);
            } else {
                auto it = name_to_camera_i.find( std::string( &strings[initial_camera_name_i] ) );
                fsc_assert( it != name_to_camera_i.end(), "could not find initial camera: " + std::string( &strings[initial_camera_name_i] ) );
                hdr->initial_camera_i = it->second;
            }
            return true;
        }

        uint field_i;
        if ( !parse_string_i( field_i, fsc, fsc_end ) ) goto error;
        const char * field = &strings[field_i];
        if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

        if ( strcmp( field, "version" ) == 0 ) {
            real version;
            if ( !parse_real( version, fsc, fsc_end, true ) ) goto error;
            fsc_assert( version == 2, fsc_file + " version != 2" );

        } else if ( strcmp( field, "camera_speed" ) == 0 ) {
            if ( !parse_real( hdr->animation_speed, fsc, fsc_end, true ) ) goto error;

        } else if ( strcmp( field, "lighting_scale" ) == 0 ) {
            if ( !parse_real( hdr->lighting_scale, fsc, fsc_end, true ) ) goto error;

        } else if ( strcmp( field, "active_camera" ) == 0 ) {
            if ( !parse_string_i( initial_camera_name_i, fsc, fsc_end ) ) goto error;

        } else if ( strcmp( field, "ambient_intensity" ) == 0 ) {
            if ( !parse_real3( hdr->ambient_intensity, fsc, fsc_end, true ) ) goto error;

        } else if ( strcmp( field, "user_defined" ) == 0 ) {
            if ( !expect_char( '{', fsc, fsc_end, true ) ) goto error;
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == '}' ) {
                    if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                    break;
                } 

                if ( !parse_string_i( field_i, fsc, fsc_end ) ) goto error;
                field = &strings[field_i];
                if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                if ( strcmp( field, "sky_box" ) == 0 ) {
                    uint name_i;
                    Texture * texture;
                    if ( !parse_string_i( name_i, fsc, fsc_end ) ) goto error;
                    if ( !load_tex( &strings[name_i], dir_name, texture ) ) goto error;
                    hdr->sky_box_tex_i = texture - textures;
                
                } else if ( strcmp( field, "sky_radius" ) == 0 ) {
                    if ( !parse_real( hdr->sky_radius, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "env_map_intensity_scale" ) == 0 ) {
                    if ( !parse_real( hdr->env_map_intensity_scale, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "opacity_scale" ) == 0 ) {
                    if ( !parse_real( hdr->opacity_scale, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "shadow_caster_count" ) == 0 ) {
                    if ( !parse_real( hdr->shadow_caster_count, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "background_color" ) == 0 ) {
                    if ( !parse_real3( hdr->background_color, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "direct_V" ) == 0 ) {
                    if ( !parse_uint( hdr->direct_V, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "direct_M" ) == 0 ) {
                    if ( !parse_uint( hdr->direct_M, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "direct_N" ) == 0 ) {
                    if ( !parse_uint( hdr->direct_N, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "tone_white" ) == 0 ) {
                    if ( !parse_real( hdr->tone_white, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "tone_key" ) == 0 ) {
                    if ( !parse_real( hdr->tone_key, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "tone_avg_luminance" ) == 0 ) {
                    if ( !parse_real( hdr->tone_avg_luminance, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "force_tone_none" ) == 0 ) {
                    if ( !parse_bool( hdr->force_tone_none, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "force_tone_avg_luminance" ) == 0 ) {
                    if ( !parse_bool( hdr->force_tone_avg_luminance, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "tex_specular_is_orm" ) == 0 ) {
                    if ( !parse_bool( hdr->tex_specular_is_orm, fsc, fsc_end ) ) goto error;

                } else {
                    fsc_assert( 0, "unexpected user_defined field '" + std::string(field) + "' in " + fsc_file );
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else if ( strcmp( field, "models" ) == 0 ) {
            if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
            uint model_file_name_i = uint(-1);
            std::string model_file_name;
            uint model_name_i = uint(-1);
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ']' ) {
                    if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                    break;
                }
                if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                for( ;; ) 
                {
                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == '}' ) {
                        if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                        break;
                    }

                    uint model_field_i;
                    if ( !parse_string_i( model_field_i, fsc, fsc_end ) ) goto error;
                    const char * model_field = &strings[model_field_i];
                    if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                    if ( strcmp( model_field, "file" ) == 0 ) {
                        if ( !parse_string_i( model_file_name_i, fsc, fsc_end ) ) goto error;
                        model_file_name = &strings[model_file_name_i];
                        if ( model_ignored_polys.find( model_file_name ) == model_ignored_polys.end() ) {
                            model_ignored_polys[model_file_name] = new std::vector<uint>;
                        }

                    } else if ( strcmp( model_field, "name" ) == 0 ) {
                        if ( !parse_string_i( model_name_i, fsc, fsc_end ) ) goto error;

                    } else if ( strcmp( model_field, "instances" ) == 0 ) {
                        if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
                        for( ;; ) 
                        {
                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ']' ) {
                                if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                                break;
                            }
                            if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                            // Instance
                            perhaps_realloc<Instance>( instances, hdr->inst_cnt, max->inst_cnt, 1 );
                            Instance * instance = &instances[hdr->inst_cnt++];
                            instance->kind = INSTANCE_KIND::MODEL;
                            instance->name_i = uint(-1);
                            instance->model_name_i = model_name_i;
                            instance->model_file_name_i = model_file_name_i;

                            // Matrix - start with identity, then compute the transpose and inverse matrices
                            perhaps_realloc<Matrix>( matrixes, hdr->matrix_cnt, max->matrix_cnt, 3 );
                            instance->matrix_i           = hdr->matrix_cnt + 0;
                            instance->matrix_inv_i       = hdr->matrix_cnt + 1;
                            instance->matrix_inv_trans_i = hdr->matrix_cnt + 2;
                            Matrix * matrix              = &matrixes[hdr->matrix_cnt++];
                            Matrix * matrix_inv          = &matrixes[hdr->matrix_cnt++];
                            Matrix * matrix_inv_trans    = &matrixes[hdr->matrix_cnt++];
                            for( uint i = 0; i < 4; i++ ) 
                            {
                                for( uint j = 0; j < 4; j++ )
                                {
                                    matrix->m[i][j] = (i == j) ? 1.0 : 0.0;
                                }
                            }

                            real3 translation = real3( 0.0, 0.0, 0.0 );
                            real3 rotation    = real3( 0.0, 0.0, 0.0 );
                            real3 scaling     = real3( 1.0, 1.0, 1.0 );
                            for( ;; ) 
                            {
                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == '}' ) {
                                    if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                                    break;
                                }

                                uint instance_field_i;
                                if ( !parse_string_i( instance_field_i, fsc, fsc_end ) ) goto error;
                                const char * instance_field = &strings[instance_field_i];
                                if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                                if ( strcmp( instance_field, "name" ) == 0 ) {
                                    if ( !parse_string_i( instance->name_i, fsc, fsc_end ) ) goto error;

                                } else if ( strcmp( instance_field, "translation" ) == 0 ) {
                                    if ( !parse_real3( translation, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( instance_field, "rotation" ) == 0 ) {
                                    if ( !parse_real3( rotation, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( instance_field, "scaling" ) == 0 ) {
                                    if ( !parse_real3( scaling, fsc, fsc_end, true ) ) goto error;
                                    matrix->scale( scaling );

                                } else {
                                    fsc_assert( 0, "unexpected instance field '" + std::string(instance_field) + "' in " + fsc_file );
                                }

                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                            }

                            // apply in reverse-TRS order
                            // rotation angles represent yaw (y-axis), pitch (x-axis), roll (z-axis), and order
                            // doesn't matter there
                            matrix->scale( scaling );
                            matrix->rotate_xz( rotation.c[0] * M_PI/180.0 );
                            matrix->rotate_yz( rotation.c[1] * M_PI/180.0 );
                            matrix->rotate_xy( rotation.c[2] * M_PI/180.0 );
                            matrix->translate( translation );

                            // pre-compute inv and inv_trans matrixes
                            matrix->invert( *matrix_inv );
                            matrix_inv->transpose( *matrix_inv_trans );

                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                        }

                    } else if ( strcmp( model_field, "ignored_polys" ) == 0 ) {
                        if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
                        fsc_assert( model_file_name_i != uint(-1), "model file name must come before ignored_polys[] in " + fsc_file );
                        std::vector<uint> * ignored_polys = model_ignored_polys[model_file_name];
                        for( ;; ) 
                        {
                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ']' ) {
                                if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                                break;
                            }

                            _int poly_i;
                            if ( !parse_int( poly_i, fsc, fsc_end ) ) goto error;
                            ignored_polys->push_back( poly_i );

                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                        }

                    } else {
                        fsc_assert( 0, "unexpected model field '" + std::string(model_field) + "' in " + fsc_file );
                    }

                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else if ( strcmp( field, "lights" ) == 0 ) {
            if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ']' ) {
                    if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                    break;
                }
                if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                // Light
                perhaps_realloc<Light>( lights, hdr->light_cnt, max->light_cnt, 1 );
                uint light_i = hdr->light_cnt;
                Light * light = &lights[hdr->light_cnt++];
                light->name_i = uint(-1);
                light->file_name_i = uint(-1);
                light->kind = LIGHT_KIND::DIRECTIONAL;
                for( uint i = 0; i < 3; i++ )
                {
                    light->intensity.c[i] = 0.0;
                    light->direction.c[i] = 0.0;
                    light->position.c[i]  = 0.0;
                    light->extent.c[i]    = 0.0;
                }
                light->opening_angle = 0.0;
                light->penumbra_angle = 0.0;
                light->diffuse_sample_cnt = 0;
                light->specular_sample_cnt = 0;

                for( ;; ) 
                {
                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == '}' ) {
                        if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                        break;
                    }

                    uint light_field_i;
                    if ( !parse_string_i( light_field_i, fsc, fsc_end ) ) goto error;
                    const char * light_field = &strings[light_field_i];
                    if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                    if ( strcmp( light_field, "name" ) == 0 ) {
                        if ( !parse_string_i( light->name_i, fsc, fsc_end ) ) goto error;
                        std::string name = &strings[light->name_i];
                        name_to_light_i[name] = light_i;

                    } else if ( strcmp( light_field, "type" ) == 0 ) {
                        uint type_name_i;
                        if ( !parse_string_i( type_name_i, fsc, fsc_end ) ) goto error;
                        const char * type_name = &strings[type_name_i];
                        if ( strcmp( type_name, "dir_light" ) == 0 ) {
                            light->kind = LIGHT_KIND::DIRECTIONAL;
                        } else if ( strcmp( type_name, "point_light" ) == 0 ) {
                            light->kind = LIGHT_KIND::POINT;
                        } else if ( strcmp( type_name, "spherical_light" ) == 0 ) {
                            light->kind = LIGHT_KIND::SPHERICAL;
                        } else if ( strcmp( type_name, "rectangular_Light" ) == 0 ) {
                            light->kind = LIGHT_KIND::RECTANGULAR;
                        } else {
                            fsc_assert( 0, "unexpected light type: " + std::string(type_name) );
                        }

                    } else if ( strcmp( light_field, "intensity" ) == 0 ) {
                        if ( !parse_real3( light->intensity, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "direction" ) == 0 ) {
                        if ( !parse_real3( light->direction, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "pos" ) == 0 ) {
                        if ( !parse_real3( light->position, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "radius" ) == 0 ) {
                        real radius;
                        if ( !parse_real( radius, fsc, fsc_end ) ) goto error;
                        light->extent.c[0] = radius;
                        light->extent.c[1] = radius;
                        light->extent.c[2] = radius;

                    } else if ( strcmp( light_field, "extent" ) == 0 ) {
                        if ( !parse_real3( light->extent, fsc, fsc_end ) ) goto error;

                    } else if ( strcmp( light_field, "opening_angle" ) == 0 ) {
                        if ( !parse_real( light->opening_angle, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "penumbra_angle" ) == 0 ) {
                        if ( !parse_real( light->penumbra_angle, fsc, fsc_end, true ) ) goto error;

                    } else {
                        fsc_assert( 0, "unexpected light field '" + std::string(light_field) + "' in " + fsc_file );
                    }

                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else if ( strcmp( field, "light_probes" ) == 0 ) {
            if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ']' ) {
                    if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                    break;
                }
                if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                // Light
                perhaps_realloc<Light>( lights, hdr->light_cnt, max->light_cnt, 1 );
                uint light_i = hdr->light_cnt;
                Light * light = &lights[hdr->light_cnt++];
                light->name_i = uint(-1);
                light->file_name_i = uint(-1);
                light->kind = LIGHT_KIND::PROBE;
                for( uint i = 0; i < 3; i++ )
                {
                    light->intensity.c[i] = 0.0;
                    light->direction.c[i] = 0.0;
                    light->position.c[i]  = 0.0;
                    light->extent.c[i]  = 0.0;
                }
                light->opening_angle = 0.0;
                light->penumbra_angle = 0.0;
                light->diffuse_sample_cnt = 0;
                light->specular_sample_cnt = 0;

                for( ;; ) 
                {
                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == '}' ) {
                        if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                        break;
                    }

                    uint light_field_i;
                    if ( !parse_string_i( light_field_i, fsc, fsc_end ) ) goto error;
                    const char * light_field = &strings[light_field_i];
                    if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                    if ( strcmp( light_field, "name" ) == 0 ) {
                        if ( !parse_string_i( light->name_i, fsc, fsc_end ) ) goto error;
                        std::string name = &strings[light->name_i];
                        name_to_light_i[name] = light_i;

                    } else if ( strcmp( light_field, "file" ) == 0 ) {
                        if ( !parse_string_i( light->file_name_i, fsc, fsc_end ) ) goto error;

                    } else if ( strcmp( light_field, "intensity" ) == 0 ) {
                        if ( !parse_real3( light->intensity, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "direction" ) == 0 ) {
                        if ( !parse_real3( light->direction, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "pos" ) == 0 ) {
                        if ( !parse_real3( light->position, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "opening_angle" ) == 0 ) {
                        if ( !parse_real( light->opening_angle, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "penumbra_angle" ) == 0 ) {
                        if ( !parse_real( light->penumbra_angle, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( light_field, "diff_samples" ) == 0 ) {
                        if ( !parse_uint( light->diffuse_sample_cnt, fsc, fsc_end ) ) goto error;

                    } else if ( strcmp( light_field, "spec_samples" ) == 0 ) {
                        if ( !parse_uint( light->specular_sample_cnt, fsc, fsc_end ) ) goto error;

                    } else {
                        fsc_assert( 0, "unexpected light probe field '" + std::string(light_field) + "' in " + fsc_file );
                    }

                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else if ( strcmp( field, "cameras" ) == 0 ) {
            if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ']' ) {
                    if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                    break;
                }
                if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                // Camera
                perhaps_realloc<Camera>( cameras, hdr->camera_cnt, max->camera_cnt, 1 );
                uint camera_i = hdr->camera_cnt;
                Camera * camera = &cameras[hdr->camera_cnt++];
                camera->name_i = uint(-1);
                for( uint i = 0; i < 3; i++ )
                {
                    camera->lookfrom.c[i] = 0.0;
                    camera->lookat.c[i] = 0.0;
                    camera->vup.c[i]  = 0.0;
                }
                camera->aperture = 0.0;
                camera->focus_dist = 21.0;
                camera->near = 0.1;
                camera->far = 10000.0;
                camera->aspect = 1.0;

                for( ;; ) 
                {
                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == '}' ) {
                        if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                        break;
                    }

                    uint camera_field_i;
                    if ( !parse_string_i( camera_field_i, fsc, fsc_end ) ) goto error;
                    const char * camera_field = &strings[camera_field_i];
                    if ( !expect_char( ':', fsc, fsc_end, true) ) goto error;

                    if ( strcmp( camera_field, "name" ) == 0 ) {
                        if ( !parse_string_i( camera->name_i, fsc, fsc_end ) ) goto error;
                        std::string name = &strings[camera->name_i];
                        name_to_camera_i[name] = camera_i;

                    } else if ( strcmp( camera_field, "pos" ) == 0 ) {
                        if ( !parse_real3( camera->lookfrom, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "target" ) == 0 ) {
                        if ( !parse_real3( camera->lookat, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "up" ) == 0 ) {
                        if ( !parse_real3( camera->vup, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "aperture" ) == 0 ) {
                        if ( !parse_real( camera->aperture, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "focal_length" ) == 0 ) {
                        if ( !parse_real( camera->focus_dist, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "depth_range" ) == 0 ) {
                        real2 near_far;
                        if ( !parse_real2( near_far, fsc, fsc_end, true ) ) goto error;
                        camera->near = near_far.c[0];
                        camera->far  = near_far.c[1];

                    } else if ( strcmp( camera_field, "vfov" ) == 0 ) {
                        if ( !parse_real( camera->vfov, fsc, fsc_end, true ) ) goto error;

                    } else if ( strcmp( camera_field, "aspect_ratio" ) == 0 ) {
                        if ( !parse_real( camera->aspect, fsc, fsc_end, true ) ) goto error;

                    } else {
                        fsc_assert( 0, "unexpected camera field '" + std::string(camera_field) + "' in " + fsc_file );
                    }

                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else if ( strcmp( field, "paths" ) == 0 ) {
            if ( !expect_char( '[', fsc, fsc_end, true ) ) goto error;
            for( ;; )
            {
                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ']' ) {
                    if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                    break;
                }
                if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                // Animation
                perhaps_realloc<Animation>( animations, hdr->animation_cnt, max->animation_cnt, 1 );
                uint animation_i = hdr->animation_cnt;
                Animation * animation = &animations[hdr->animation_cnt++];
                animation->name_i = uint(-1);
                animation->start_camera_i = uint(-1);
                animation->start_frame_i = hdr->frame_cnt;
                animation->frame_cnt = 0;
                animation->is_looped = false;

                for( ;; ) 
                {
                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == '}' ) {
                        if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                        break;
                    }

                    uint animation_field_i;
                    if ( !parse_string_i( animation_field_i, fsc, fsc_end ) ) goto error;
                    const char * animation_field = &strings[animation_field_i];
                    if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                    if ( strcmp( animation_field, "name" ) == 0 ) {
                        if ( !parse_string_i( animation->name_i, fsc, fsc_end ) ) goto error;
                        std::string name = &strings[animation->name_i];
                        name_to_animation_i[name] = animation_i;

                    } else if ( strcmp( animation_field, "loop" ) == 0 ) {
                        if ( !parse_bool( animation->is_looped, fsc, fsc_end ) ) goto error;

                    } else if ( strcmp( animation_field, "attached_objects" ) == 0 ) {
                        skip_whitespace( fsc, fsc_end );
                        if ( !expect_char( '[', fsc, fsc_end ) ) goto error;
                        for( ;; ) 
                        {
                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ']' ) {
                                if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                                break;
                            }
                            if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                            for( ;; ) 
                            {
                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == '}' ) {
                                    if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                                    break;
                                }

                                uint attached_field_i;
                                if ( !parse_string_i( attached_field_i, fsc, fsc_end ) ) goto error;
                                const char * attached_field = &strings[attached_field_i];
                                if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                                if ( strcmp( attached_field, "type" ) == 0 ) {
                                    uint type_i;
                                    if ( !parse_string_i( type_i, fsc, fsc_end ) ) goto error;
                                    const char * type = &strings[type_i];
                                    fsc_assert( strcmp( type, "camera" ) == 0, "attached_objects type must be camera for now" );

                                } else if ( strcmp( attached_field, "name" ) == 0 ) {
                                    if ( !parse_string_i( animation->start_camera_i, fsc, fsc_end ) ) goto error;
                                
                                } else {
                                    fsc_assert( 0, "unexpected attached_objects field '" + std::string(attached_field) + "' in " + fsc_file );
                                }

                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                            }

                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                        }

                    } else if ( strcmp( animation_field, "frames" ) == 0 ) {
                        skip_whitespace( fsc, fsc_end );
                        if ( !expect_char( '[', fsc, fsc_end ) ) goto error;
                        for( ;; ) 
                        {
                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ']' ) {
                                if ( !expect_char( ']', fsc, fsc_end ) ) goto error;
                                break;
                            }
                            if ( !expect_char( '{', fsc, fsc_end ) ) goto error;

                            // Frame
                            perhaps_realloc<Frame>( frames, hdr->frame_cnt, max->frame_cnt, 1 );
                            Frame * frame = &frames[hdr->frame_cnt++];
                            frame->time = 0.0;
                            frame->camera_i = hdr->camera_cnt++;

                            // Camera
                            perhaps_realloc<Camera>( cameras, hdr->camera_cnt, max->camera_cnt, 1 );
                            Camera * camera = &cameras[frame->camera_i];
                            std::string camera_name = &strings[animation->name_i];
                            camera_name += "_" + std::to_string( animation->frame_cnt - 1 );
                            uint cn_len = camera_name.length();
                            perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, cn_len+1 );
                            char * to_cn = &strings[hdr->char_cnt];
                            hdr->char_cnt += cn_len + 1;
                            memcpy( to_cn, camera_name.c_str(), cn_len+1 );
                            for( uint i = 0; i < 3; i++ )
                            {
                                camera->lookfrom.c[i] = 0.0;
                                camera->lookat.c[i] = 0.0;
                                camera->vup.c[i]  = 0.0;
                            }
                            camera->aperture = -1.0;
                            camera->focus_dist = -1.0;
                            camera->near = -1.0;
                            camera->far = -1.0;
                            camera->vfov = 40.0;
                            camera->aspect = -1.0;

                            for( ;; ) 
                            {
                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == '}' ) {
                                    if ( !expect_char( '}', fsc, fsc_end ) ) goto error;
                                    break;
                                }

                                uint frame_field_i;
                                if ( !parse_string_i( frame_field_i, fsc, fsc_end ) ) goto error;
                                const char * frame_field = &strings[frame_field_i];
                                if ( !expect_char( ':', fsc, fsc_end, true ) ) goto error;

                                if ( strcmp( frame_field, "time" ) == 0 ) {
                                    if ( !parse_real( frame->time, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "pos" ) == 0 ) {
                                    if ( !parse_real3( camera->lookfrom, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "target" ) == 0 ) {
                                    if ( !parse_real3( camera->lookat, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "up" ) == 0 ) {
                                    if ( !parse_real3( camera->vup, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "aperture" ) == 0 ) {
                                    if ( !parse_real( camera->aperture, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "focal_length" ) == 0 ) {
                                    if ( !parse_real( camera->focus_dist, fsc, fsc_end, true ) ) goto error;

                                } else if ( strcmp( frame_field, "depth_range" ) == 0 ) {
                                    real2 near_far;
                                    if ( !parse_real2( near_far, fsc, fsc_end, true ) ) goto error;
                                    camera->near = near_far.c[0];
                                    camera->far  = near_far.c[1];

                                } else if ( strcmp( frame_field, "aspect_ratio" ) == 0 ) {
                                    if ( !parse_real( camera->aspect, fsc, fsc_end, true ) ) goto error;

                                } else {
                                    fsc_assert( 0, "unexpected frame field '" + std::string(frame_field) + "' in " + fsc_file );
                                }

                                skip_whitespace( fsc, fsc_end );
                                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                            }

                            skip_whitespace( fsc, fsc_end );
                            fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                            if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                        }

                    } else {
                        fsc_assert( 0, "unexpected frame field '" + std::string(animation_field) + "' in " + fsc_file );
                    }

                    skip_whitespace( fsc, fsc_end );
                    fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                    if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
                }

                skip_whitespace( fsc, fsc_end );
                fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
                if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
            }

        } else {
            fsc_assert( 0, "unexpected fscene field '" + std::string(field) + "' in " + fsc_file );
        }

        skip_whitespace( fsc, fsc_end );
        fsc_assert( fsc != fsc_end, "unexpected end of " + fsc_file );
        if ( *fsc == ',' && !expect_char( ',', fsc, fsc_end ) ) goto error;
    }

    error:
        error_msg += " (at line " + std::to_string( line_num ) + " of " + fsc_file + ")";
        die_assert( false, "aborting..." );
        return false;
}

bool Model::load_obj( std::string obj_file, std::string dir_name )
{
    (void)dir_name;

    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    line_num = 1;
    if ( !file_read( obj_file, obj_start, obj_end ) ) return false;
    obj = obj_start;

    //------------------------------------------------------------
    // Parse .obj file contents
    //------------------------------------------------------------
    const char *mtllib = nullptr;
    const char *obj_name = nullptr;
    const char *name;
    std::string mtl_name;
    Material *  material = nullptr;
    uint        mtl_i = uint(-1);
    Object *    object = nullptr;
    Polygon *   polygon = nullptr;
    Vertex *    vertex = nullptr;

    for( ;; ) 
    {
        skip_whitespace( obj, obj_end );
        if ( obj == obj_end ) return true;

        obj_cmd_t cmd;
        if ( !parse_obj_cmd( cmd ) ) break;
        dprint( "obj_cmd=" + std::to_string( cmd ) );

        switch( cmd )
        {
            case CMD_O:
                perhaps_realloc<Object>( objects, hdr->obj_cnt, max->obj_cnt, 1 );
                object = &objects[ hdr->obj_cnt++ ];

                if ( !parse_name( obj_name, obj, obj_end ) ) goto error;
                object->name_i = obj_name - strings;
                object->poly_cnt = 0;
                object->poly_i = hdr->poly_cnt;
                break;
                
            case CMD_G:
                if ( !parse_name( name, obj, obj_end ) ) goto error;
                break;
                
            case CMD_S:
                if ( !parse_name( name, obj, obj_end ) ) goto error;
                break;
                
            case CMD_V:
                perhaps_realloc<real3>( positions, hdr->pos_cnt, max->pos_cnt, 1 );
                if ( !parse_real3( positions[ hdr->pos_cnt++ ], obj, obj_end ) ) goto error;
                skip_to_eol( obj, obj_end );
                break;
                
            case CMD_VN:
                perhaps_realloc<real3>( normals, hdr->norm_cnt, max->norm_cnt, 1 );
                if ( !parse_real3( normals[ hdr->norm_cnt++ ], obj, obj_end ) ) goto error;
                skip_to_eol( obj, obj_end );
                break;
                
            case CMD_VT:
                perhaps_realloc<real2>( texcoords, hdr->texcoord_cnt, max->texcoord_cnt, 1 );
                if ( !parse_real2( texcoords[ hdr->texcoord_cnt++ ], obj, obj_end ) ) goto error;
                skip_to_eol( obj, obj_end );
                break;
                
            case CMD_F:
            {
                perhaps_realloc<Polygon>( polygons, hdr->poly_cnt, max->poly_cnt, 1 );
                uint poly_i = hdr->poly_cnt++;
                polygon = &polygons[poly_i];
                polygon->mtl_i = mtl_i;
                polygon->vtx_cnt = 0;
                polygon->vtx_i = hdr->vtx_cnt;

                while( !eol( obj, obj_end ) ) 
                {
                    polygon->vtx_cnt++;
                    perhaps_realloc<Vertex>( vertexes, hdr->vtx_cnt, max->vtx_cnt, 1 );
                    vertex = &vertexes[ hdr->vtx_cnt++ ];

                    _int v_i;
                    if ( !parse_int( v_i, obj, obj_end ) )   goto error;
                    dprint( "v_i=" + std::to_string( v_i ) );
                    vertex->v_i = (v_i >= 0)  ? v_i : (hdr->pos_cnt + v_i);

                    if ( obj != obj_end && *obj == '/' ) {
                        if ( !expect_char( '/', obj, obj_end ) ) goto error;
                        _int vt_i;
                        if ( obj != obj_end && *obj == '/' ) {
                            if ( !expect_char( '/', obj, obj_end ) ) goto error;
                            vertex->vt_i = uint(-1);
                        } else {
                            if ( !parse_int( vt_i, obj, obj_end ) ) goto error;
                            dprint( "vt_i=" + std::to_string( vt_i ) );
                            vertex->vt_i = (vt_i >= 0) ? vt_i : ((int)hdr->texcoord_cnt + vt_i);
                        }

                        if ( obj != obj_end && *obj == '/' ) {
                            if ( !expect_char( '/', obj, obj_end ) ) goto error;
                            _int vn_i;
                            if ( !parse_int( vn_i, obj, obj_end ) )  goto error;
                            dprint( "vn_i=" + std::to_string( vn_i ) );
                            vertex->vn_i = (vn_i >= 0) ? vn_i : (hdr->norm_cnt + vn_i);
                        } else {
                            vertex->vn_i = uint(-1);
                        }
                    } else {
                        vertex->vt_i = uint(-1);
                        vertex->vn_i = uint(-1);
                    }
                }
                obj_assert( polygon->vtx_cnt != 0, ".obj f command has no vertices" );
                if ( object != nullptr ) object->poly_cnt++;
                // split 4+ sided polygons into triangles
                uint first_poly_i = hdr->poly_cnt-1;
                if ( polygon->vtx_cnt > 3 ) {
                    uint extra_triangle_cnt = polygon->vtx_cnt - 3;
                    if ( 0 ) {
                        Vertex * pvertexes = &vertexes[polygon->vtx_i];
                        real3 p0( positions[pvertexes[0].v_i].c[0], positions[pvertexes[0].v_i].c[1], positions[pvertexes[0].v_i].c[2] );
                        real3 p1( positions[pvertexes[1].v_i].c[0], positions[pvertexes[1].v_i].c[1], positions[pvertexes[1].v_i].c[2] );
                        real3 p2( positions[pvertexes[2].v_i].c[0], positions[pvertexes[2].v_i].c[1], positions[pvertexes[2].v_i].c[2] );
                        real3 p3( positions[pvertexes[3].v_i].c[0], positions[pvertexes[3].v_i].c[1], positions[pvertexes[3].v_i].c[2] );
                        std::cout << "splitting poly_i=" << first_poly_i << " " << *polygon <<
                                     " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << " p3=" << p3 << "\n";
                    }
                    uint poly_vtx_i = polygon->vtx_i;

                    // save all polygon vertexes here (to make this easier)
                    obj_assert( polygon->vtx_cnt < 128, "too many vertices" );
                    Vertex saved_vertexes[128];
                    for( uint v = 0; v < polygon->vtx_cnt; v++ ) saved_vertexes[v] = vertexes[poly_vtx_i+v];

                    // allocate more triangles
                    perhaps_realloc<Polygon>( polygons, hdr->poly_cnt, max->poly_cnt, extra_triangle_cnt );
                    polygon = polygons + first_poly_i;  // could have changed
                    hdr->poly_cnt += extra_triangle_cnt;

                    // allocate 2*extra_triangle_cnt more vertexes
                    perhaps_realloc<Vertex>( vertexes, hdr->vtx_cnt, max->vtx_cnt, 2*extra_triangle_cnt );
                    hdr->vtx_cnt += 2*extra_triangle_cnt;

                    // fill in triangles as a triangle fan
                    for( uint t = 0; t < (1+extra_triangle_cnt); t++ )
                    {
                        Polygon * triangle = polygon + t;
                        triangle->mtl_i    = polygon->mtl_i;
                        triangle->vtx_i    = poly_vtx_i + 3*t;
                        triangle->vtx_cnt  = 3;
                        vertexes[triangle->vtx_i+0] = saved_vertexes[0];    // common point
                        vertexes[triangle->vtx_i+1] = saved_vertexes[t+1];
                        vertexes[triangle->vtx_i+2] = saved_vertexes[t+2];
                    }
                }
                // precompute surface normal and area for each triangle (ignore for points and lines)
                uint triangle_cnt = (polygon->vtx_cnt < 3) ? 0 : (hdr->poly_cnt - first_poly_i);
                for( uint t = 0; t < triangle_cnt; t++ )
                {
                    Polygon * triangle = polygon + t;
                    Vertex * pvertexes = &vertexes[triangle->vtx_i];
                    real3 p0( positions[pvertexes[0].v_i].c[0], positions[pvertexes[0].v_i].c[1], positions[pvertexes[0].v_i].c[2] );
                    real3 p1( positions[pvertexes[1].v_i].c[0], positions[pvertexes[1].v_i].c[1], positions[pvertexes[1].v_i].c[2] );
                    real3 p2( positions[pvertexes[2].v_i].c[0], positions[pvertexes[2].v_i].c[1], positions[pvertexes[2].v_i].c[2] );
                    triangle->normal = (p1 - p0).cross( p2 - p0 );
                    real len = triangle->normal.length();
                    triangle->area = len / 2;
                    triangle->normal = triangle->normal / len;
                    if ( false && triangle_cnt > 1 ) std::cout << "    poly_i=" << first_poly_i+t << " " << *triangle << ": p0=" << p0 << " p1=" << p1 << " p2=" << p2 << "\n";

                    if ( mtl_i != uint(-1) && materials[mtl_i].is_emissive() ) {
                        perhaps_realloc<uint>( emissive_polygons, hdr->emissive_poly_cnt, max->emissive_poly_cnt, 1 );
                        emissive_polygons[hdr->emissive_poly_cnt++] = first_poly_i + t;
                    }
                }
                break;
            }

            case CMD_MTLLIB:
                if ( !parse_name( mtllib, obj, obj_end ) ) goto error;
                if ( !load_mtl( mtllib, dir_name ) ) goto error;
                break;
                
            case CMD_USEMTL:
                obj_assert( mtllib != nullptr, "no mtllib defined for object " + std::string( obj_name ) );
                if ( !eol( obj, obj_end ) ) {
                    if ( !parse_name( name, obj, obj_end ) ) goto error;
                    mtl_name = std::string( name );
                    obj_assert( name_to_mtl_i.find( mtl_name ) != name_to_mtl_i.end(), "unknown material: " + std::string( mtl_name ) );
                    mtl_i = name_to_mtl_i[mtl_name];
                    material = &materials[mtl_i];
                } else {
                    mtl_i = uint(-1);
                    material = nullptr;
                }
                break;

            default:
                break;
        }
    }
    error:
        error_msg += " (at line " + std::to_string( line_num ) + " of " + obj_file + ")";
        die_assert( false, "aborting..." );
        return false;
}

bool Model::load_mtl( std::string mtl_file, std::string dir_name, bool replacing )
{
    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    if ( dir_name != std::string( "" ) ) mtl_file = dir_name + "/" + mtl_file;
    if ( !file_read( mtl_file, mtl_start, mtl_end ) ) return false;
    mtl = mtl_start;

    //------------------------------------------------------------
    // Parse .mtl file contents
    //------------------------------------------------------------
    uint orig_mtl_cnt = hdr->mtl_cnt;
    if ( replacing ) hdr->mtl_cnt = 0;

    const char * mtl_name = nullptr;
    const char * tex_name = nullptr;

    uint        mtl_i;
    Material *  material = nullptr;
    uint        tex_i;
    Texture *   texture  = nullptr;
    std::string option_name;
    real        bump_multiplier;
    real3       tex_offset;
    real3       tex_scale;

    for( ;; ) 
    {
        skip_whitespace( mtl, mtl_end );
        if ( mtl == mtl_end ) break;

        mtl_cmd_t cmd;
        if ( !parse_mtl_cmd( cmd, mtl, mtl_end ) ) return false;
        dprint( "  mtl_cmd: " + std::to_string( cmd ) );

        switch( cmd )
        {
            case CMD_NEWMTL:
                if ( !replacing ) {
                    if ( !parse_name( mtl_name, mtl, mtl_end ) ) return false;
                } else {
                    skip_to_eol( mtl, mtl_end );
                }
                if ( !replacing && name_to_mtl_i.find( mtl_name ) != name_to_mtl_i.end() ) {
                    mtl_i = name_to_mtl_i[ mtl_name ];         
                    material = &materials[mtl_i];
                    dprint( "  found " + std::string( mtl_name ) );
                } else {
                    if ( !replacing ) perhaps_realloc<Material>( materials, hdr->mtl_cnt, max->mtl_cnt, 1 );
                    rtn_assert( !replacing || hdr->mtl_cnt < orig_mtl_cnt, "more replacement materials than originally created" );
                    mtl_i = hdr->mtl_cnt++;
                    material = &materials[mtl_i];
                    if ( !replacing ) {
                        material->name_i = mtl_name - strings;
                        for( uint j = 0; j < 3; j++ )
                        {
                            material->Ka.c[j] = 1.0;
                            material->Kd.c[j] = 1.0;
                            material->Ke.c[j] = 0.0;
                            material->Ks.c[j] = 0.0;
                            material->Tf.c[j] = 1.0;
                        }
                        material->Tr = 0.0;
                        material->Ns = 0.0;
                        material->Ni = 1.0;
                        material->d  = 1.0;
                        material->illum = 2;
                        material->map_Ka_i = uint(-1);
                        material->map_Kd_i = uint(-1);
                        material->map_Ke_i = uint(-1);
                        material->map_Ks_i = uint(-1);
                        material->map_Ns_i = uint(-1);
                        material->map_d_i  = uint(-1);
                        material->map_Bump_i = uint(-1);
                        material->map_refl_i = uint(-1);
                        name_to_mtl_i[ mtl_name ] = mtl_i;
                        dprint( "  added " + std::string( mtl_name ) );
                    }
                }
                break;

            case CMD_KA:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real3( material->Ka, mtl, mtl_end ) ) return false;
                break;

            case CMD_KD:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real3( material->Kd, mtl, mtl_end ) ) return false;
                break;

            case CMD_KE:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real3( material->Ke, mtl, mtl_end ) ) return false;
                break;

            case CMD_KS:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real3( material->Ks, mtl, mtl_end ) ) return false;
                break;

            case CMD_TF:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real3( material->Tf, mtl, mtl_end ) ) return false;
                break;

            case CMD_TR:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real( material->Tr, mtl, mtl_end ) ) return false;
                material->d = 1.0 - material->Tr;
                break;

            case CMD_NS:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real( material->Ns, mtl, mtl_end ) ) return false;
                break;

            case CMD_NI:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real( material->Ni, mtl, mtl_end ) ) return false;
                break;

            case CMD_D:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real( material->d, mtl, mtl_end ) ) return false;
                break;

            case CMD_ILLUM:
                rtn_assert( material != nullptr, "no material defined" );
                if ( !parse_real( material->illum, mtl, mtl_end ) ) return false;
                break;

            case CMD_MAP_KA:
            case CMD_MAP_KD:
            case CMD_MAP_KE:
            case CMD_MAP_KS:
            case CMD_MAP_NS:
            case CMD_MAP_D:
            case CMD_MAP_BUMP:
            case CMD_MAP_REFL:
                rtn_assert( material != nullptr, "no material defined" );
                
                if ( !replacing ) {
                    bump_multiplier = 1.0;
                    tex_offset = real3( 0, 0, 0 );
                    tex_scale  = real3( 1, 1, 1 );
                    while( parse_option_name( option_name, mtl, mtl_end ) )
                    {
                        if ( option_name == std::string("-bm") ) {
                            if ( !parse_real( bump_multiplier, mtl, mtl_end ) ) return false;
                        } else if ( option_name == std::string("-o") ) {
                            if ( !parse_real3( tex_offset, mtl, mtl_end ) ) return false;
                        } else if ( option_name == std::string("-s") ) {
                            if ( !parse_real3( tex_scale, mtl, mtl_end ) ) return false;
                        } else {
                            rtn_assert( 0, "unknown texture map option: " + option_name );
                        }
                    }
                    if ( !parse_name( tex_name, mtl, mtl_end ) ) return false;
                    if ( name_to_tex_i.find( tex_name ) != name_to_tex_i.end() ) {
                        // already loaded it
                        //
                        tex_i = name_to_tex_i[tex_name];
                        texture = &textures[tex_i];
                    } else {
                        // need to load it
                        //
                        if ( !load_tex( tex_name, dir_name, texture ) ) return false;
                    }
                    texture->bump_multiplier = bump_multiplier;
                    texture->offset          = tex_offset;
                    texture->scale           = tex_scale;

                    tex_i = texture - textures;

                    switch( cmd ) 
                    {
                        case CMD_MAP_KA:        material->map_Ka_i   = tex_i; break;
                        case CMD_MAP_KD:        material->map_Kd_i   = tex_i; break;
                        case CMD_MAP_KE:        material->map_Ke_i   = tex_i; break;
                        case CMD_MAP_KS:        material->map_Ks_i   = tex_i; break;
                        case CMD_MAP_NS:        material->map_Ns_i   = tex_i; break;
                        case CMD_MAP_D:         material->map_d_i    = tex_i; break;
                        case CMD_MAP_BUMP:      material->map_Bump_i = tex_i; break;
                        case CMD_MAP_REFL:      material->map_refl_i = tex_i; break;
                        default:                                              break; // should not happen
                    }
                } else {
                    skip_to_eol( mtl, mtl_end );
                }
                break;

            default:
                rtn_assert( 0, "unknown .mtl command" );
        }

        rtn_assert( eol( mtl, mtl_end ), "not at eol in .mtl file: " + surrounding_lines( mtl, mtl_end ) );
    }

    return true;
}

uint Model::debug_tex_i = 1; // uint(-1);

bool Model::load_tex( const char * tex_name, std::string dir_name, Texture *& texture, const unsigned char * data, uint w, uint h, uint nchan )
{
    //---------------------------------------------
    // Allocate Texture structure
    //---------------------------------------------
    perhaps_realloc<Texture>( textures, hdr->tex_cnt, max->tex_cnt, 1 );
    uint tex_i = hdr->tex_cnt++;
    texture = &textures[tex_i];
    memset( texture, 0, sizeof( Texture ) );
    texture->name_i = tex_name - strings;
    name_to_tex_i[ tex_name ] = tex_i;

    //---------------------------------------------
    // Change any '\' to '/'
    //---------------------------------------------
    std::string tex_name_s = std::string( tex_name );
    if ( dir_name != "" ) tex_name_s = dir_name + "/" + tex_name_s;
    char * file_name = strdup( tex_name_s.c_str() );
    size_t len = strlen( file_name );
    for( size_t i = 0; i < len; i++ ) 
    {
        if ( file_name[i] == '\\' ) file_name[i] = '/';
    }

    //---------------------------------------------
    // Dissect the path so we can replace the extension.
    //---------------------------------------------
    bool from_file = data == nullptr;
    std::string base_name;
    std::string ext_name;
    dissect_path( file_name, dir_name, base_name, ext_name );

    //---------------------------------------------
    // There are three possible regimes:
    // 1) read non-ASTC file for miplevel 0, generate other miplevels, keep uncompressed
    // 2) read non-ASTC file for miplevel 0, generate other miplevels, 
    //    write out uncompressed file, compress to ASTC using ARM's astcenc program,
    //    read back in .astc file.  In this case, we need to keep uncompressed
    //    texels in a separate buffer.
    // 3) read ASTC file for each miplevel.  All must exist.
    //---------------------------------------------
    std::string astc_name = dir_name + "/" + base_name + ".astc";
    bool        astc_file_exists     = file_exists( astc_name );
    bool        in_uncompressed_mode = hdr->texture_compression == TEXTURE_COMPRESSION::NONE;
    bool        in_astc_gen_mode     = !astc_file_exists && hdr->texture_compression == TEXTURE_COMPRESSION::ASTC;
    bool        in_astc_read_mode    =  astc_file_exists && hdr->texture_compression == TEXTURE_COMPRESSION::ASTC;

    texture->compression = hdr->texture_compression;
    if ( texture->compression == TEXTURE_COMPRESSION::ASTC ) {
        //---------------------------------------------
        // Align to 16B boundary to match ASTC alignment
        // Should not need to reallocate due to >16B granularity of allocations.
        //---------------------------------------------
        hdr->texel_cnt += ((hdr->texel_cnt % 16) == 0) ? 0 : (16 - (hdr->texel_cnt % 16)); 
    }
    texture->texel_i = hdr->texel_cnt;  // always

    uint width = 0;
    uint height = 0;
    uint byte_width = 0;
    uint byte_cnt = 0;
    for( uint mipmap_level = 0; ; mipmap_level++ )
    {
        if ( !in_astc_read_mode ) {
            if ( mipmap_level == 0 ) {
                //---------------------------------------------
                // READ NON-ASTC FILE FOR MIP 0
                //---------------------------------------------
                if ( data == nullptr ) {
                    int ww, hh, nnchan;
                    data = stbi_load( file_name, &ww, &hh, &nnchan, 0 );
                    rtn_assert( data != nullptr, "unable to read in texture file " + std::string( file_name ) + ", error reason: " + std::string(stbi_failure_reason()) );
                    w = ww;
                    h = hh;
                    nchan = nnchan;
                } else {
                    assert( w > 0 && h > 0 && nchan > 0 );
                }
                texture->width  = w;
                texture->height = h;
                texture->nchan  = nchan;
                width      = texture->width;
                height     = texture->height;
                byte_width = width * texture->nchan;
                byte_cnt   = byte_width * height;

            } else if ( mipmap_level != 0 && !in_astc_read_mode ) {
                //---------------------------------------------
                // GENERATE UNCOMPRESSED TEXELS FOR THIS MIP LEVEL
                //
                // To keep this simple, we create a new buffer, 
                // then copy it back to texels[] later if we're in_uncompressed_mode.
                //---------------------------------------------
                uint to_width  = width  >> 1;
                uint to_height = height >> 1;
                if ( to_width  == 0 ) to_width  = 1;
                if ( to_height == 0 ) to_height = 1;

                uint   to_byte_width = to_width * texture->nchan;
                uint   to_byte_cnt   = to_byte_width * to_height;
                unsigned char * to_data = aligned_alloc<unsigned char>( to_byte_cnt );
                unsigned char * trgb = to_data;
                for( uint ti = 0; ti < to_height; ti++ )
                {
                    for( uint tj = 0; tj < to_width; tj++, trgb += texture->nchan )
                    {
                        uint sum[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
                        uint cnt = 0;
                        for( uint fi = 0; fi < 2; fi++ )
                        {
                            uint fii = ti*2 + fi;
                            if ( fii >= height ) continue;
                            for( uint fj = 0; fj < 2; fj++ )
                            {
                                uint fjj = tj*2 + fj;
                                if ( fjj >= width ) continue;
                                const unsigned char * frgb = data + fii*byte_width + fjj*texture->nchan;
                                for( uint c = 0; c < texture->nchan; c++ )
                                {
                                    sum[c] += frgb[c];
                                }
                                cnt++;
                            }
                        }

                        for( uint c = 0; c < texture->nchan; c++ )
                        {
                            trgb[c] = (real(sum[c]) + 0.5) / real(cnt);
                        }
                    }
                }
                delete data;

                height          = to_height;
                width           = to_width;
                byte_width      = to_byte_width;
                byte_cnt        = to_byte_cnt;
                data            = to_data;
            }

        } else if ( mipmap_level != 0 ) {
            //---------------------------------------------
            // Update width and height.
            //---------------------------------------------
            if ( width  != 1 ) width  >>= 1;
            if ( height != 1 ) height >>= 1;
        }

        if ( in_uncompressed_mode ) {
            //---------------------------------------------
            // SAVE UNCOMPRESSED TEXELS
            //
            // For consistency, keep the data buffer for mipmap generation in the next pass.
            //---------------------------------------------
            perhaps_realloc<unsigned char>( texels, hdr->texel_cnt, max->texel_cnt, byte_cnt );
            memcpy( texels + hdr->texel_cnt, data, byte_cnt );
            hdr->texel_cnt += byte_cnt;

        } else if ( in_astc_gen_mode ) {
            //---------------------------------------------
            // WRITE OUT UNCOMPRESSED TEXELS AND CREATE ASTC FILE
            //
            // Write out temporary .BMP file, then 
            // use ARM's astcenc program to compress the texels.
            // It must be on the current PATH.
            //---------------------------------------------
            int mypid = getpid();
            std::string bmp_name = dir_name + "/" + base_name + "." + std::to_string(mipmap_level) + 
                                   "." + std::to_string(mypid) + ".bmp";
            if ( !stbi_write_bmp( bmp_name.c_str(), width, height, texture->nchan, data ) ) {
                die_assert( false, "could not write out temporary BMP file " + bmp_name );    
            }

            std::string suff = (mipmap_level == 0) ? "" : ("." + std::to_string(mipmap_level));
            std::string astc_name = dir_name + "/" + base_name + suff + ".astc";
            uint        dim_x;
            uint        dim_y;
            Texture::astc_blk_dims_get( width, height, dim_x, dim_y ); // chooses some default
            cmd( "astcenc -c '" + bmp_name + "' '" + astc_name + "' " + std::to_string(dim_x) + "x" + std::to_string(dim_y) + " -exhaustive > /dev/null" );
            cmd( "rm -f '" + bmp_name + "'" );
            die_assert( file_exists( astc_name ), "weird, ASTC file " + astc_name + " did not get created by astcenc" );
        } 

        if ( in_astc_gen_mode || in_astc_read_mode ) {
            //---------------------------------------------
            // READ ASTC FILE 
            //
            // It should exist at this point.
            //---------------------------------------------
            std::string suff = (mipmap_level == 0) ? "" : ("." + std::to_string(mipmap_level));
            astc_name = dir_name + "/" + base_name + suff + ".astc";
            die_assert( file_exists( astc_name ), "ASTC file " + astc_name + " does not exist for miplevel " + 
                                                  std::to_string(mipmap_level) + "; delete them all and try again" );
            char * start;
            char * end;
            if ( !file_read( astc_name, start, end ) ) {
                die_assert( false, "unable to read ASTC file for miplevel " + std::to_string(mipmap_level) + " even though file exists: " + astc_name );
            }
            unsigned char * astc_data = reinterpret_cast<unsigned char *>( start );
            uint astc_byte_cnt = end - start;

            //---------------------------------------------
            // We assume ARM .astc file format where there is a 16B header block 
            // with the following format.  Little endian.
            //---------------------------------------------
            die_assert( astc_byte_cnt >= 32, "ASTC file " + astc_name + " does not have at least 32 bytes (header + one block)" );
            die_assert( (astc_byte_cnt % 16) == 0, "ASTC file " + astc_name + " does not have a multiple of 16 bytes" );

            const Texture::ASTC_Header * astc_hdr = reinterpret_cast<const Texture::ASTC_Header *>( astc_data );
            die_assert( astc_hdr->magic_is_good(), "ASTC file " + astc_name + " header magic number is not expected value" );
            die_assert( astc_hdr->blockdim_z == 1, "ASTC file " + astc_name + " uses 3D textures which are currently unsuported" );
            uint w = (astc_hdr->xsize[0] << 0) | (astc_hdr->xsize[1] << 8) | (astc_hdr->xsize[2] << 16);
            uint h = (astc_hdr->ysize[0] << 0) | (astc_hdr->ysize[1] << 8) | (astc_hdr->ysize[2] << 16);
            uint d = (astc_hdr->zsize[0] << 0) | (astc_hdr->zsize[1] << 8) | (astc_hdr->zsize[2] << 16);
            die_assert( d == 1, "ASTC file " + astc_name + " uses 3D textures which are currently unsupported" );
            if ( mipmap_level == 0 ) {
                texture->width  = w;
                texture->height = h;
                texture->nchan  = 4;
                width  = w;
                height = h;
            } else {
                die_assert( w == width,  "unexpected width in ASTC file " + astc_name );
                die_assert( h == height, "unexpected height in ASTC file " + astc_name );
            }
            uint64_t x_blk_cnt = (w + astc_hdr->blockdim_x - 1) / astc_hdr->blockdim_x;
            uint64_t y_blk_cnt = (h + astc_hdr->blockdim_y - 1) / astc_hdr->blockdim_y;
            uint64_t file_blk_cnt = astc_byte_cnt/16 - 1;
            die_assert( file_blk_cnt == (x_blk_cnt * y_blk_cnt), "ASTC file " + astc_name + " does not have the expected number of blocks" );

            //---------------------------------------------
            // Use the same header to keep aligned to 16B.
            // We memcpy() the entire file - as is - to our texels[] space,
            // aligning on 16B boundary.
            //---------------------------------------------
            perhaps_realloc<unsigned char>( texels, hdr->texel_cnt, max->texel_cnt, astc_byte_cnt );
            memcpy( texels + hdr->texel_cnt, astc_data, astc_byte_cnt );
            hdr->texel_cnt += astc_byte_cnt;
            delete astc_data;
        }

        //---------------------------------------------
        // Stop if not generating mip levels > 0 or we've just done the 1x1 level.
        //---------------------------------------------
        if ( hdr->mipmap_filter == MIPMAP_FILTER::NONE || (width == 1 && height == 1) ) break;
    } 
    if ( from_file && data != nullptr ) delete data;
    return true;
}

bool Model::load_gph( std::string gph_file, std::string dir_name, std::string base_name, std::string ext_name )
{
    (void)dir_name;
    (void)base_name;
    (void)ext_name;

    //------------------------------------------------------------
    // Map in .gph file
    //------------------------------------------------------------
    line_num = 1;
    if ( !file_read( gph_file, gph_start, gph_end ) ) return false;
    gph = gph_start;

    uint block_i = gph_node_alloc( GRAPH_NODE_KIND::BLOCK );
    hdr->graph_root_i = block_i;

    //------------------------------------------------------------
    // Parse .gph file contents
    //------------------------------------------------------------
    uint id_i;
    uint op_i;
    uint expr_i;
    for( ;; ) 
    {
        skip_whitespace( gph, gph_end );
        if ( gph == gph_end ) return true;

        //---------------------------------------------
        // Parse assign
        //---------------------------------------------
        if ( !parse_gph_id( id_i, gph, gph_end ) ) return false;
        if ( !parse_gph_op( op_i, gph, gph_end ) ) return false;
        die_assert( strcmp( &strings[graph_nodes[op_i].u.s_i], "=" ) == 0, "bad assignment operator" );
        if ( !parse_gph_expr( expr_i, gph, gph_end ) ) return false;

        //---------------------------------------------
        // Add ASSIGN node
        //---------------------------------------------
        uint assign_i = gph_node_alloc( GRAPH_NODE_KIND::ASSIGN, block_i );
        graph_nodes[assign_i].u.child_first_i = id_i;
        graph_nodes[id_i].sibling_i = expr_i;
        graph_nodes[id_i].parent_i = assign_i;
        graph_nodes[expr_i].parent_i = assign_i;
    }
    return true;
}

uint Model::gph_node_alloc( Model::GRAPH_NODE_KIND kind, uint parent_i )
{
    //mdout << "new node: kind=" << kind << "\n";
    perhaps_realloc<Graph_Node>( graph_nodes, hdr->graph_node_cnt, max->graph_node_cnt, 1 );
    uint node_i = hdr->graph_node_cnt++;
    Graph_Node * node = &graph_nodes[node_i];
    node->kind = kind;
    node->u.child_first_i = uint(-1);
    node->sibling_i = uint(-1);

    // optionally add to tail of parent's child list
    node->parent_i = parent_i;
    if ( parent_i != uint(-1) ) {
        Graph_Node * parent = &graph_nodes[parent_i];
        for( uint * sibling_i_ptr = &parent->u.child_first_i; ; sibling_i_ptr = &graph_nodes[*sibling_i_ptr].sibling_i )
        {
            if ( *sibling_i_ptr == uint(-1) ) {
                *sibling_i_ptr = node_i;
                break;
            }
        }
    }
    return node_i;
}

bool Model::load_nvdb( std::string nvdb_file, std::string dir_name, std::string base_name, std::string ext_name )
{
    (void)dir_name;
    (void)base_name;
    (void)ext_name;

    //------------------------------------------------------------
    // Map in .nvdb file
    //------------------------------------------------------------
    line_num = 1;
    if ( !file_read( nvdb_file, nvdb_start, nvdb_end ) ) return false;
    nvdb = nvdb_start;

    perhaps_realloc<Volume>( volumes, hdr->volume_cnt, max->volume_cnt, 1 );
    uint volume_i = hdr->volume_cnt++;
    Volume * volume = &volumes[volume_i];
    volume->grid_cnt = 0;
    volume->grid_i = hdr->volume_grid_cnt; // first

    //------------------------------------------------------------
    // Parse file.
    //------------------------------------------------------------
    while( nvdb != nvdb_end ) 
    {
        //------------------------------------------------------------
        // Parse next segment header.
        //------------------------------------------------------------
        enum class Codec : uint16_t { NONE = 0, ZIP = 1, BLOSC = 2, END = 3 };

        struct SegmentHeader  // in the file
        {
            uint64_t    magic;
            uint16_t    major;
            uint16_t    minor;
            uint16_t    grid_cnt;
            Codec       codec;
        };
        die_assert( sizeof(SegmentHeader) == 16, "unexpected SegmentHeader size" );
        die_assert( size_t(nvdb_end-nvdb) >= sizeof(SegmentHeader), "can't read SegmentHeader - not enough bytes in .nvdb file" );
        SegmentHeader header;
        memcpy( &header, nvdb, sizeof(SegmentHeader) );
        nvdb += sizeof(SegmentHeader);
        die_assert( header.magic == 0x304244566f6e614eUL, "bad magic number in .nvdb segment header for grid_i=" + str(volume->grid_cnt));
        die_assert( header.major == 19, "bad major number in .nvdb segment header" );
        volume->grid_cnt += header.grid_cnt;

        const VolumeGrid * grid0 = &volume_grids[volume->grid_i];
        for( uint i = 0; i < header.grid_cnt; i++ )
        {
            //------------------------------------------------------------
            // Parse next GridMetaData.
            //------------------------------------------------------------
            struct AABB64
            {
                real64          _min[3];
                real64          _max[3];
            };

            struct AABBI32
            {
                int32_t         _min[3];
                int32_t         _max[3];
            };

            struct GridMetaData   // in the file
            {
                uint64_t        grid_size; 
                uint64_t        file_size; 
                uint64_t        name_key; 
                uint64_t        voxel_cnt;
                VolumeVoxelType grid_type;
                VolumeGridClass grid_class;
                uint32_t        name_size;      // includes \0
                uint32_t        node_cnt[4];
                AABB64          world_box;
                AABBI32         index_box;
                double          voxel_size;     // can't use this; need to look at Map
            };
            die_assert( sizeof(GridMetaData) == 144, "unexpected GridMetaData size" );
            die_assert( size_t(nvdb_end-nvdb) >= sizeof(GridMetaData), "can't read GridMetaData - not enough bytes in .nvdb file" );
            GridMetaData meta;
            memcpy( &meta, nvdb, sizeof(GridMetaData) );
            nvdb += sizeof(GridMetaData);

            die_assert( meta.name_size > 0, "name_size=0 in GridMetaData in .nvdb file" );
            std::string name( nvdb );
            nvdb += meta.name_size;
            uint name_i = make_string( name );

            //------------------------------------------------------------
            // Allocate new VolumeGrid and fill in header information.
            //------------------------------------------------------------
            perhaps_realloc<VolumeGrid>( volume_grids, hdr->volume_grid_cnt, max->volume_grid_cnt, 1 );
            uint grid_i = hdr->volume_grid_cnt++;
            VolumeGrid * grid = &volume_grids[grid_i];
            grid->name_i = name_i;
            grid->voxel_cnt = meta.voxel_cnt;
            grid->voxel_type = meta.grid_type;
            grid->voxel_class = (name == "d" || name == "dens"     || name == "density")      ? VolumeVoxelClass::DENSITY :
                                (name == "t" || name == "temp"     || name == "temperature")  ? VolumeVoxelClass::TEMPERATURE : 
                                (name == "r" || name == "rgb"      || name == "RGB")          ? VolumeVoxelClass::RGB : 
                                (name == "e" || name == "emissive" || name == "emissive_rgb") ? VolumeVoxelClass::EMISSIVE_RGB : 
                                (name == "n" || name == "norm"     || name == "normal")       ? VolumeVoxelClass::NORMAL : 
                                (name == "i" || name == "ior"      || name == "IOR")          ? VolumeVoxelClass::IOR : 
                                (name == "f" || name == "f0"       || name == "F0")           ? VolumeVoxelClass::F0 : 
                                (name == "g" || name == "G")                                  ? VolumeVoxelClass::G : 
                                (name == "a" || name == "A"        || name == "attenuation")  ? VolumeVoxelClass::ATTENUATION : 
                                (name == "m" || name == "M"        || name == "MFPL" || name == "mfpl") ? VolumeVoxelClass::MFPL : 
                                                                                                VolumeVoxelClass::UNKNOWN;
            grid->grid_class = meta.grid_class;
            for( uint j = 0; j < 4; j++ ) 
            { 
                grid->node_cnt[j] = meta.node_cnt[j];
                if ( j < 3 ) {
                    grid->world_box._min.c[j] = meta.world_box._min[j];
                    grid->world_box._max.c[j] = meta.world_box._max[j];
                    grid->index_box._min[j]   = meta.index_box._min[j];
                    grid->index_box._max[j]   = meta.index_box._max[j];
                    if ( i > 0 ) {
                        die_assert( grid->world_box._min.c[j] == grid0->world_box._min.c[j], "all grids in a volume must have same world box" );
                        die_assert( grid->world_box._max.c[j] == grid0->world_box._max.c[j], "all grids in a volume must have same world box" );
                        die_assert( grid->index_box._min[j]   == grid0->index_box._min[j],   "all grids in a volume must have same index box" );
                        die_assert( grid->index_box._max[j]   == grid0->index_box._max[j],   "all grids in a volume must have same index box" );
                    }
                }
            }

            //------------------------------------------------------------
            // Copy voxel data which is variable size.
            // We always align this on an 8B boundary in the voxels[] array.
            // Each voxels[] element is one byte even though each voxel in our
            // grid will typically consume multiple bytes.
            //
            // Note that the voxel data is actually an elaborate tree structure
            // that matches the NanoVDB format. We don't interpret the structure here.
            // That occurs in the real_value() et al. methods that take voxel [x,y,z] 
            // and walk the tree to get the required information.  The tree walk 
            // is non-trivial.
            //
            // In other words, the voxel data is not a simple 3D array of floats. :-)
            //------------------------------------------------------------
            die_assert( header.codec == Codec::NONE, ".nvdb files must be uncompressed for now" );
            die_assert( meta.file_size == meta.grid_size, "uncompressed grids must have equal file_size and grid_size" );
            uint64_t alloc_voxel_cnt = ((meta.grid_size+7)/8) * 8;
            perhaps_realloc<unsigned char>( voxels, hdr->voxel_cnt, max->voxel_cnt, alloc_voxel_cnt );
            grid->voxel_i = hdr->voxel_cnt;
            hdr->voxel_cnt += alloc_voxel_cnt;

            memcpy( &voxels[grid->voxel_i], nvdb, meta.grid_size );
            nvdb += meta.file_size;

            //------------------------------------------------------------
            // Obtain the world_voxel_size and world_translate from the NVDB MapData
            // within the GridData.
            //------------------------------------------------------------
            struct MapData { 
                float  mMatF[9]; 
                float  mInvMatF[9];
                float  mVecF[3];
                float  mTaperF;
                double mMatD[9];
                double mInvMatD[9];
                double mVecD[3];
                double mTaperD;
            };

            static const int MaxNameSize = 256;
            struct GridData { 
                uint64_t    mMagic;             // 8 byte magic to validate it is valid grid data.
                                                // REFERENCEMAGIC = 0x4244566f6e614eL;
                char        mGridName[MaxNameSize];
                real64      mBBoxMin[3];        // min corner of AABB
                real64      mBBoxMax[3];        // max corner of AABB
                MapData     mMap;               // affine transformation between index and world space in both single and double precision
                double      mUniformScale;      // size of a voxel in world units
                VolumeGridClass mGridClass;     // 4 bytes
                VolumeVoxelType mVoxelType;     // 4 bytes
                uint32_t    mBlindDataCount;    // count of GridBlindMetaData structures that follow this grid (after the gridname).
                uint32_t    padding[2];         // seems to be needed to match up with NanoVDB
            };

            const GridData * ngrid = reinterpret_cast< const GridData * >( &voxels[grid->voxel_i] );

            grid->world_voxel_size = real3d( ngrid->mMap.mMatD[0], ngrid->mMap.mMatD[4], ngrid->mMap.mMatD[8] );
            grid->world_voxel_size_inv = real3d(1,1,1) / grid->world_voxel_size;
            grid->world_translate  = real3d( ngrid->mMap.mVecD[0], ngrid->mMap.mVecD[1], ngrid->mMap.mVecD[2] );
            std::cout << "grids[" << grid_i << "] voxel_class=" << grid->voxel_class << " world_box=" << grid->world_box << 
                         " world_voxel_size=" << grid->world_voxel_size << " world_voxel_size_inv=" << grid->world_voxel_size_inv << 
                         " world_translate=" << grid->world_translate << "\n";

            if ( i != 0 ) {
                die_assert( grid->world_voxel_size == grid0->world_voxel_size, "all grids in same volume must have the same world_voxel_size" );
                die_assert( grid->world_translate  == grid0->world_translate,  "all grids in same volume must have the same world_translate" );
            }
        }
    }
    return true;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// TEXTURE QUERIES
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real4 Model::Texture::texel_read( const Model * model, uint mip_level, uint64 ui, uint64 vi, uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const
{
    switch( compression )
    {
        case TEXTURE_COMPRESSION::NONE: return texel_read_uncompressed( model, mip_level, ui, vi, vaddr, byte_cnt, do_srgb_to_linear );
        case TEXTURE_COMPRESSION::ASTC: return texel_read_astc(         model, mip_level, ui, vi, vaddr, byte_cnt, do_srgb_to_linear );
        default:
        {
            std::cout << "ERROR: unexpected texture compression: " << int(compression) << "\n";
            die_assert( false, "aborting..." );
            return real4( 0, 0, 0, 0 );
        }
    }
}

inline Model::real4 Model::Texture::texel_read_uncompressed( const Model * model, uint mip_level, uint64 ui, uint64 vi, uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const
{
    const unsigned char * mdata = &model->texels[texel_i];
    uint mw = width;
    uint mh = height;
    if ( model->hdr->mipmap_filter != MIPMAP_FILTER::NONE ) {
        // find the beginning of the proper mip texture 
        for( _int imip_level = mip_level; imip_level > 0 && !(mw == 1 && mh == 1); imip_level-- ) 
        {
            mdata += nchan * mw * mh;
            if ( mw != 1 ) mw >>= 1;
            if ( mh != 1 ) mh >>= 1;
        }
    }

    uint64_t toffset = nchan*ui + nchan*mw*vi;

    if (vaddr != nullptr)    *vaddr = reinterpret_cast<uint64_t>( &mdata[toffset] );
    if (byte_cnt != nullptr) *byte_cnt = 3;

    uint ru = mdata[toffset+0];
    uint gu = (nchan  > 2) ? mdata[toffset+1] : ru;
    uint bu = (nchan  > 2) ? mdata[toffset+2] : ru;
    uint au = (nchan  > 3) ? mdata[toffset+3] : (nchan == 2) ? mdata[toffset+1] : 255;

    real r = do_srgb_to_linear ? srgb8_to_linear_gamma( ru ) : (real(ru) / 255.0);
    real g = do_srgb_to_linear ? srgb8_to_linear_gamma( gu ) : (real(gu) / 255.0);
    real b = do_srgb_to_linear ? srgb8_to_linear_gamma( bu ) : (real(bu) / 255.0);
    real a = real(au) / 255.0;  // always stored linear
    return real4( r, g, b, a );
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// COLOR SPACE CONVERSION UTILITIES
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real Model::srgb_to_linear_gamma( Model::real srgb )
{
    if ( srgb <= 0.04045 ) {
        return srgb / 12.92;
    } else {
        return std::pow( (srgb + 0.055)/1.055, 2.4 );
    }
}

inline Model::real Model::srgb8_to_linear_gamma( uint8_t srgb )
{
    static bool did_lut_init = false;
    static real lut[256];
    if ( !did_lut_init ) {
        for( uint i = 0; i < 256; i++ )
        {
            lut[i] = srgb_to_linear_gamma( real(i) / 255.0 );
        } 
        did_lut_init = true;
    }
    return lut[srgb];
}

Model::real Model::linear_to_srgb_gamma( Model::real linear )
{
    if ( linear < 0.0031308 ) {
        return linear * 12.92;
    } else {
        return 1.055 * std::pow( linear, 1.0/2.4 ) - 0.055;
    }
}

Model::real Model::linear8_to_srgb_gamma( uint8_t linear )
{
    static bool did_lut_init = false;
    static real lut[256];
    if ( !did_lut_init ) {
        for( uint i = 0; i < 256; i++ )
        {
            lut[i] = linear_to_srgb_gamma( real(i) / 255.0 );
        } 
        did_lut_init = true;
    }
    return lut[linear];
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// VECTOR (2D, 3D, 4D) AND MATRIX (3D)
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real Model::real4::dot( const Model::real4 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1] + c[2] * v2.c[2] + c[3] * v2.c[3];
}

inline Model::real Model::real4::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3] ); 
}

inline Model::real Model::real4::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3];
}

inline Model::real Model::real4::distance( const Model::real4& v2 ) const 
{ 
    auto diff = *this - v2;
    return diff.length();
}

inline Model::real Model::real4::distance_sqr( const Model::real4& v2 ) const 
{ 
    auto diff = *this - v2;
    return diff.length_sqr();
}

inline Model::real4& Model::real4::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real4 Model::real4::normalized( void ) const
{
    return *this / length();
}

inline void Model::real4::clamp( real min, real max ) 
{
    if ( c[0] < min ) c[0] = 0;
    if ( c[1] < min ) c[1] = 0;
    if ( c[2] < min ) c[2] = 0;
    if ( c[3] < min ) c[3] = 0;
    if ( c[0] > max ) c[0] = 1;
    if ( c[1] > max ) c[1] = 1;
    if ( c[2] > max ) c[2] = 1;
    if ( c[3] > max ) c[3] = 1;
}

inline Model::real4 Model::real4::clamped() const 
{
    Model::real4 v = *this;
    v.clamp();
    return v;
}

inline Model::real4 Model::real4::mulby2() const 
{
    return Model::real4( ::mulby2(c[0]), ::mulby2(c[1]), ::mulby2(c[2]), ::mulby2(c[3]) );
}

inline Model::real4 Model::real4::mulby4() const 
{
    return Model::real4( ::mulby4(c[0]), ::mulby4(c[1]), ::mulby4(c[2]), ::mulby4(c[3]) );
}

inline Model::real4 Model::real4::divby2() const 
{
    return Model::real4( ::divby2(c[0]), ::divby2(c[1]), ::divby2(c[2]), ::divby2(c[3]) );
}

inline Model::real4 Model::real4::divby4() const 
{
    return Model::real4( ::divby4(c[0]), ::divby4(c[1]), ::divby4(c[2]), ::divby4(c[3]) );
}

inline bool Model::real4::operator == ( const Model::real4 &v2 ) const
{
    return c[0] == v2.c[0] && c[1] == v2.c[1] && c[2] == v2.c[2] && c[3] == v2.c[3];
}

inline Model::real4& Model::real4::operator += ( const Model::real4 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    c[2] += v2.c[2];
    c[3] += v2.c[3];
    return *this;
}

inline Model::real4& Model::real4::operator -= ( const Model::real4 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    c[2] -= v2.c[2];
    c[3] -= v2.c[3];
    return *this;
}

inline Model::real4& Model::real4::operator *= ( const Model::real4 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    c[2] *= v2.c[2];
    c[3] *= v2.c[3];
    return *this;
}

inline Model::real4& Model::real4::operator *= ( const Model::real s )
{
    c[0] *= s;
    c[1] *= s;
    c[2] *= s;
    c[3] *= s;
    return *this;
}

inline Model::real4& Model::real4::operator /= ( const Model::real4 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    c[2] /= v2.c[2];
    c[3] /= v2.c[3];
    return *this;
}

inline Model::real4& Model::real4::operator /= ( const Model::real s )
{
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
    c[3] /= s;
    return *this;
}

inline Model::real64 Model::real3d::dot( const Model::real3d &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1] + c[2] * v2.c[2];
}

inline Model::real3d Model::real3d::cross( const Model::real3d &v2 ) const
{
    return real3d( (c[1]*v2.c[2]   - c[2]*v2.c[1]),
                  (-(c[0]*v2.c[2] - c[2]*v2.c[0])),
                  (c[0]*v2.c[1]   - c[1]*v2.c[0]) );
}

inline Model::real64 Model::real3d::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] ); 
}

inline Model::real64 Model::real3d::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
}

inline Model::real64 Model::real3d::distance( const Model::real3d& v2 ) const 
{ 
    auto diff = *this - v2;
    return diff.length();
}

inline Model::real64 Model::real3d::distance_sqr( const Model::real3d& v2 ) const 
{ 
    auto diff = *this - v2;
    return diff.length_sqr();
}

inline Model::real3d& Model::real3d::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real3d Model::real3d::normalized( void ) const
{
    return *this / length();
}

inline bool Model::real3d::operator == ( const Model::real3d &v2 ) const
{
    return c[0] == v2.c[0] && c[1] == v2.c[1] && c[2] == v2.c[2];
}

inline void Model::real3d::clamp( real64 min, real64 max ) 
{
    if ( c[0] < min ) c[0] = 0;
    if ( c[1] < min ) c[1] = 0;
    if ( c[2] < min ) c[2] = 0;
    if ( c[0] > max ) c[0] = 1;
    if ( c[1] > max ) c[1] = 1;
    if ( c[2] > max ) c[2] = 1;
}

inline Model::real3d Model::real3d::clamped( void ) const 
{
    Model::real3d v = *this;
    v.clamp();
    return v;
}

inline Model::real3d Model::real3d::mulby2( void ) const 
{
    return Model::real3d( ::mulby2(c[0]), ::mulby2(c[1]), ::mulby2(c[2]) );
}

inline Model::real3d Model::real3d::mulby4( void ) const 
{
    return Model::real3d( ::mulby4(c[0]), ::mulby4(c[1]), ::mulby4(c[2]) );
}

inline Model::real3d Model::real3d::divby2( void ) const 
{
    return Model::real3d( ::divby2(c[0]), ::divby2(c[1]), ::divby2(c[2]) );
}

inline Model::real3d Model::real3d::divby4( void ) const 
{
    return Model::real3d( ::divby4(c[0]), ::divby4(c[1]), ::divby4(c[2]) );
}

inline Model::real3d& Model::real3d::operator += ( const Model::real3d &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    c[2] += v2.c[2];
    return *this;
}

inline Model::real3d& Model::real3d::operator -= ( const Model::real3d &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    c[2] -= v2.c[2];
    return *this;
}

inline Model::real3d& Model::real3d::operator *= ( const Model::real3d &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    c[2] *= v2.c[2];
    return *this;
}

inline Model::real3d& Model::real3d::operator *= ( const Model::real64 s )
{
    c[0] *= s;
    c[1] *= s;
    c[2] *= s;
    return *this;
}

inline Model::real3d& Model::real3d::operator /= ( const Model::real3d &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    c[2] /= v2.c[2];
    return *this;
}

inline Model::real3d& Model::real3d::operator /= ( const Model::real64 s )
{
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
    return *this;
}

inline Model::real Model::real3::dot( const Model::real3 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1] + c[2] * v2.c[2];
}

inline Model::real3 Model::real3::cross( const Model::real3 &v2 ) const
{
    return real3( (c[1]*v2.c[2]   - c[2]*v2.c[1]),
                  (-(c[0]*v2.c[2] - c[2]*v2.c[0])),
                  (c[0]*v2.c[1]   - c[1]*v2.c[0]) );
}

inline Model::real Model::real3::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] ); 
}

inline Model::real Model::real3::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
}

inline Model::real Model::real3::distance_sqr( const Model::real3& v2 ) const 
{ 
    auto diff = *this - v2;
    return diff.length_sqr();
}

inline Model::real3& Model::real3::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real3 Model::real3::normalized( void ) const
{
    return *this / length();
}

inline void Model::real3::clamp( real min, real max ) 
{
    if ( c[0] < min ) c[0] = 0;
    if ( c[1] < min ) c[1] = 0;
    if ( c[2] < min ) c[2] = 0;
    if ( c[0] > max ) c[0] = 1;
    if ( c[1] > max ) c[1] = 1;
    if ( c[2] > max ) c[2] = 1;
}

inline Model::real3 Model::real3::clamped( void ) const 
{
    Model::real3 v = *this;
    v.clamp();
    return v;
}

inline Model::real3 Model::real3::mulby2( void ) const 
{
    return Model::real3( ::mulby2(c[0]), ::mulby2(c[1]), ::mulby2(c[2]) );
}

inline Model::real3 Model::real3::mulby4( void ) const 
{
    return Model::real3( ::mulby4(c[0]), ::mulby4(c[1]), ::mulby4(c[2]) );
}

inline Model::real3 Model::real3::divby2( void ) const 
{
    return Model::real3( ::divby2(c[0]), ::divby2(c[1]), ::divby2(c[2]) );
}

inline Model::real3 Model::real3::divby4( void ) const 
{
    return Model::real3(::divby4(c[0]), ::divby4(c[1]), ::divby4(c[2]));
}

inline bool Model::real3::operator == ( const Model::real3 &v2 ) const
{
    return c[0] == v2.c[0] && c[1] == v2.c[1] && c[2] == v2.c[2];
}

inline Model::real3& Model::real3::operator += ( const Model::real3 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    c[2] += v2.c[2];
    return *this;
}

inline Model::real3& Model::real3::operator -= ( const Model::real3 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    c[2] -= v2.c[2];
    return *this;
}

inline Model::real3& Model::real3::operator *= ( const Model::real3 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    c[2] *= v2.c[2];
    return *this;
}

inline Model::real3& Model::real3::operator *= ( const Model::real s )
{
    c[0] *= s;
    c[1] *= s;
    c[2] *= s;
    return *this;
}

inline Model::real3& Model::real3::operator /= ( const Model::real3 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    c[2] /= v2.c[2];
    return *this;
}

inline Model::real3& Model::real3::operator /= ( const Model::real s )
{
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
    return *this;
}

inline Model::real Model::real2::dot( const Model::real2 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1];
}

inline Model::real Model::real2::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] );
}

inline Model::real Model::real2::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1];
}

inline Model::real2& Model::real2::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real2 Model::real2::normalized( void ) const
{
    return *this / length();
}

inline void Model::real2::clamp( real min, real max ) 
{
    if ( c[0] < min ) c[0] = 0;
    if ( c[1] < min ) c[1] = 0;
    if ( c[0] > max ) c[0] = 1;
    if ( c[1] > max ) c[1] = 1;
}

inline Model::real2 Model::real2::clamped( void ) const 
{
    Model::real2 v = *this;
    v.clamp();
    return v;
}

inline Model::real2& Model::real2::operator += ( const Model::real2 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    return *this;
}

inline Model::real2& Model::real2::operator -= ( const Model::real2 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    return *this;
}

inline Model::real2& Model::real2::operator *= ( const Model::real2 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    return *this;
}

inline Model::real2& Model::real2::operator *= ( const Model::real s )
{
    c[0] *= s;
    c[1] *= s;
    return *this;
}

inline Model::real2& Model::real2::operator /= ( const Model::real2 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    return *this;
}

inline Model::real2& Model::real2::operator /= ( const Model::real s )
{
    c[0] /= s;
    c[1] /= s;
    return *this;
}

inline Model::Matrix Model::Matrix::make_view( const real3& lookfrom, const real3& lookat, const real3& vup )
{
    real3 forward = lookat - lookfrom;
    forward.normalize();

    real3 side = forward.cross( vup );
    side.normalize();

    real3 up = side.cross( forward );
    
    Model::Matrix M;
    M.m[0][0] = side[0];
    M.m[0][1] = side[1];
    M.m[0][2] = side[2];
    M.m[0][3] = 0.0;
    M.m[1][0] = up[0];
    M.m[1][1] = up[1];
    M.m[1][2] = up[2];
    M.m[1][3] = 0.0;
    M.m[2][0] = -forward[0];
    M.m[2][1] = -forward[1];
    M.m[2][2] = -forward[2];
    M.m[2][3] = 0.0;
    M.m[3][0] = 0.0;
    M.m[3][1] = 0.0;
    M.m[3][2] = 0.0;
    M.m[3][3] = 1.0;

    M.translate( -lookfrom );

    return M;
}

inline Model::Matrix Model::Matrix::make_frustum( real left, real right, real bottom, real top, real near, real far )
{
    real near_2x = 2.0 * near;
    real right_m_left = right - left;
    real top_m_bottom = top - bottom;
    real far_m_near = far - near;  
    Model::Matrix M;
    M.m[0][0] = near_2x / right_m_left;
    M.m[0][1] = 0.0;
    M.m[0][2] = (right + left) / right_m_left;
    M.m[0][3] = 0.0;
    M.m[1][0] = 0.0;
    M.m[1][1] = near_2x / top_m_bottom;
    M.m[1][2] = (top + bottom) / top_m_bottom;
    M.m[1][3] = 0.0;
    M.m[2][0] = 0.0;
    M.m[2][1] = 0.0;
    M.m[2][2] = -(far + near) / far_m_near;
    M.m[2][3] = -(far * near_2x) / far_m_near;
    M.m[3][0] = 0.0;
    M.m[3][1] = 0.0;
    M.m[3][2] = -1.0;
    M.m[3][3] = 0.0;
    return M;
}

inline Model::Matrix Model::Matrix::make_perspective( real vfov, real aspect, real near, real far )
{
    real bottom = -near * std::tan( vfov * PI_DIV_360 );
    real top    = -bottom;
    real left   = aspect * bottom;
    real right  = -left;
    return make_frustum( left, right, bottom, top, near, far );
}

inline bool Model::Matrix::is_identity( void ) const
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            if (m[i][j] != ((i == j) ? 1 : 0)) return false;
        }
    }
    return true;
}

inline void Model::Matrix::identity( void )
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            m[i][j] = (i == j) ? 1 : 0;
        }
    }
}

inline void Model::Matrix::translate( const real3& translation )
{
    m[0][3] += translation.c[0];
    m[1][3] += translation.c[1];
    m[2][3] += translation.c[2];
}

inline void Model::Matrix::scale( const real3& scaling )
{
    m[0][0] *= scaling.c[0];
    m[1][1] *= scaling.c[1];
    m[2][2] *= scaling.c[2];
}

void Model::Matrix::rotate_xy( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[0][0] = c;
    M2.m[0][1] = s;
    M2.m[1][0] = -s;
    M2.m[1][1] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

void Model::Matrix::rotate_xz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[0][0] = c;
    M2.m[0][2] = s;
    M2.m[2][0] = -s;
    M2.m[2][2] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

void Model::Matrix::rotate_yz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[1][1] = c;
    M2.m[1][2] = -s;
    M2.m[2][1] = s;
    M2.m[2][2] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

inline Model::Matrix Model::Matrix::operator + ( const Matrix& m ) const
{
    Matrix r;
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            r.m[i][j] = this->m[i][j] + m.m[i][j];
        }
    }
    return r;
}

inline Model::Matrix Model::Matrix::operator - ( const Matrix& m ) const
{
    Matrix r;
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            r.m[i][j] = this->m[i][j] - m.m[i][j];
        }
    }
    return r;
}

inline bool Model::Matrix::operator == ( const Matrix& m ) const
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            if ( this->m[i][j] != m.m[i][j] ) return false;
        }
    }
    return true;
}

inline void Model::Matrix::multiply( double s ) 
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            m[i][j] *= s;
        }
    }
}

Model::real4 Model::Matrix::row( uint r ) const
{
    real4 v;
    for( uint32_t c = 0; c < 4; c++ ) 
    {
        v.c[c] = m[r][c];
    }
    return v;
}

Model::real4 Model::Matrix::column( uint c ) const
{
    real4 v;
    for( uint32_t r = 0; r < 4; r++ ) 
    {
        v.c[r] = m[r][c];
    }
    return v;
}

inline void Model::Matrix::transform( const real4& v, real4& r ) const
{
    if (is_identity()) {
        r = v;
    } else {
        // order: r = *this * v
        for( uint i = 0; i < 4; i++ )
        {
            double sum = 0.0;               // use higher-precision here
            for( uint j = 0; j < 4; j++ )
            {
                double partial = m[i][j];
                partial *= v.c[j];
                sum += partial;
            }
            r.c[i] = sum;
        }
    }
}

void Model::Matrix::transform( const real3& v, real3& r, bool div_by_w ) const
{
    // order: r = *this * v
    if ( div_by_w ) {
        real4 v4 = real4( v.c[0], v.c[1], v.c[2], 1.0 );
        real4 r4;
        transform( v4, r4 );
        r.c[0] = r4.c[0];
        r.c[1] = r4.c[1];
        r.c[2] = r4.c[2];
        r /= r4.c[3];                   // w
    } else {
        for( uint i = 0; i < 3; i++ )
        {
            double sum = 0.0;               // use higher-precision here
            for( uint j = 0; j < 3; j++ )
            {
                double partial = m[i][j];
                partial *= v.c[j];
                sum += partial;
            }
            r.c[i] = sum;
        }
    }
}

void Model::Matrix::transform( const real3d& v, real3d& r, bool div_by_w ) const
{
    // order: r = *this * v
    if ( div_by_w ) {
        die( "can't transform real3d with div_by_w=true right now, can if needed" );
    } else {
        for( uint i = 0; i < 3; i++ )
        {
            double sum = 0.0;               // use higher-precision here
            for( uint j = 0; j < 3; j++ )
            {
                double partial = m[i][j];
                partial *= v.c[j];
                sum += partial;
            }
            r.c[i] = sum;
        }
    }
}

inline void Model::Matrix::transform( const Ray& r, Ray& r2 ) const
{
    real3 origin2;
    real3 direction2;
    transform( r.origin(), origin2, false );
    for( uint i = 0; i < 3; i++ ) 
    {
        origin2.c[i] += m[i][3];
    }
    transform( r.direction(), direction2, false );
    r2 = Ray( origin2, direction2, r.kind(), r.solid_angle(), r.cone_angle() );
}

inline void Model::Matrix::transform( const Ray64& r, Ray64& r2 ) const
{
    real3d origin2;
    real3d direction2;
    transform( r.origin(), origin2, false );
    for( uint i = 0; i < 3; i++ ) 
    {
        origin2.c[i] += m[i][3];
    }
    transform( r.direction(), direction2, false );
    r2 = Ray64( origin2, direction2, r.kind(), r.solid_angle(), r.cone_angle() );
}

void Model::Matrix::transform( const Matrix& M2, Matrix& M3 ) const
{
    // order: M3 = *this * M2
    for( uint r = 0; r < 4; r++ )
    {
        for( uint c = 0; c < 4; c++ )
        {
            double sum = 0.0;
            for( int k = 0; k < 4; k++ )
            {
                double partial = m[r][k];
                partial *= M2.m[k][c];
                sum += partial;
            }
            M3.m[r][c] = sum;
        }
    }
}

void Model::Matrix::transpose( Model::Matrix& mt ) const
{
    for( int i = 0; i < 4; i++ )
    {
        for( int j = 0; j < 4; j++ )
        {
            mt.m[j][i] = m[i][j];
        }
    }
}

//---------------------------------------------------------------
// Invert a 4x4 Matrix that is known to be an affine transformation.
// This is simpler and faster, but does not work for perspective matrices and other cases.
//---------------------------------------------------------------
void Model::Matrix::invert_affine( Model::Matrix& minv ) const
{
    double i00 = m[0][0];
    double i01 = m[0][1];
    double i02 = m[0][2];
    double i03 = m[0][3];
    double i10 = m[1][0];
    double i11 = m[1][1];
    double i12 = m[1][2];
    double i13 = m[1][3];
    double i20 = m[2][0];
    double i21 = m[2][1];
    double i22 = m[2][2];
    double i23 = m[2][3];
    double i30 = m[3][0];
    double i31 = m[3][1];
    double i32 = m[3][2];
    double i33 = m[3][3];

    double s0  = i00 * i11 - i10 * i01;
    double s1  = i00 * i12 - i10 * i02;
    double s2  = i00 * i13 - i10 * i03;
    double s3  = i01 * i12 - i11 * i02;
    double s4  = i01 * i13 - i11 * i03;
    double s5  = i02 * i13 - i12 * i03;

    double c5  = i22 * i33 - i32 * i23;
    double c4  = i21 * i33 - i31 * i23;
    double c3  = i21 * i32 - i31 * i22;
    double c2  = i20 * i33 - i30 * i23;
    double c1  = i20 * i32 - i30 * i22;
    double c0  = i20 * i31 - i30 * i21;

    double invdet  = 1 / ( s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0 );

    minv.m[0][0] = (i11 * c5 - i12 * c4 + i13 * c3) * invdet;
    minv.m[0][1] = (-i01 * c5 + i02 * c4 - i03 * c3) * invdet;
    minv.m[0][2] = (i31 * s5 - i32 * s4 + i33 * s3) * invdet;
    minv.m[0][3] = (-i21 * s5 + i22 * s4 - i23 * s3) * invdet;

    minv.m[1][0] = (-i10 * c5 + i12 * c2 - i13 * c1) * invdet;
    minv.m[1][1] = (i00 * c5 - i02 * c2 + i03 * c1) * invdet;
    minv.m[1][2] = (-i30 * s5 + i32 * s2 - i33 * s1) * invdet;
    minv.m[1][3] = (i20 * s5 - i22 * s2 + i23 * s1) * invdet;

    minv.m[2][0] = (i10 * c4 - i11 * c2 + i13 * c0) * invdet;
    minv.m[2][1] = (-i00 * c4 + i01 * c2 - i03 * c0) * invdet;
    minv.m[2][2] = (i30 * s4 - i31 * s2 + i33 * s0) * invdet;
    minv.m[2][3] = (-i20 * s4 + i21 * s2 - i23 * s0) * invdet;

    minv.m[3][0] = (-i10 * c3 + i11 * c1 - i12 * c0) * invdet;
    minv.m[3][1] = (i00 * c3 - i01 * c1 + i02 * c0) * invdet;
    minv.m[3][2] = (-i30 * s3 + i31 * s1 - i32 * s0) * invdet;
    minv.m[3][3] = (i20 * s3 - i21 * s1 + i22 * s0) * invdet;
}

//---------------------------------------------------------------
// This is a more complex invert() but works for all types of matrices.
//---------------------------------------------------------------
void Model::Matrix::invert( Model::Matrix& minv ) const 
{
    // Inverse = adjoint / determinant
    adjoint( minv );

    // Determinant is the dot product of the first row and the first row
    // of cofactors (i.e. the first col of the adjoint matrix)
    real4 col0 = minv.column( 0 );
    real4 row0 = minv.row( 0 );
    double det = col0.dot( row0 );
    minv.multiply( 1.0/det );
}

void Model::Matrix::adjoint( Model::Matrix& M ) const 
{
    Matrix Mt;     
    cofactor( Mt );
    Mt.transpose( M );
}

double Model::Matrix::determinant( void ) const 
{
    // Determinant is the dot product of the first row and the first row
    // of cofactors (i.e. the first col of the adjoint matrix)
    Matrix C;
    cofactor( C );
    real4 row0 = row( 0 );
    return row0.dot( row0 );
}

void Model::Matrix::cofactor( Model::Matrix& C ) const 
{
    // We'll use i to incrementally compute -1 ^ (r+c)
    int32_t i = 1;

    for( uint r = 0; r < 4; r++ ) 
    {
        for( uint c = 0; c < 4; c++ ) 
        {
            // Compute the determinant of the 3x3 submatrix
            double det = subdeterminant( r, c );
            C.m[r][c] = double(i) * det;
            i = -i;
        }
        i = -i;
    }
}

double Model::Matrix::subdeterminant( uint exclude_row, uint exclude_col ) const 
{
    // Compute non-excluded row and column indices
    uint _row[3];
    uint _col[3];

    for( uint i = 0; i < 3; i++ ) 
    {
        _row[i] = i;
        _col[i] = i;

        if( i >= exclude_row ) _row[i]++;
        if( i >= exclude_col ) _col[i]++;
    }

    // Compute the first row of cofactors 
    double cofactor00 =
      double(m[_row[1]][_col[1]]) * double(m[_row[2]][_col[2]]) -
      double(m[_row[1]][_col[2]]) * double(m[_row[2]][_col[1]]);

    double cofactor10 =
      double(m[_row[1]][_col[2]]) * double(m[_row[2]][_col[0]]) -
      double(m[_row[1]][_col[0]]) * double(m[_row[2]][_col[2]]);

    double cofactor20 =
      double(m[_row[1]][_col[0]]) * double(m[_row[2]][_col[1]]) -
      double(m[_row[1]][_col[1]]) * double(m[_row[2]][_col[0]]);

    // Product of the first row and the cofactors along the first row
    return
      double(m[_row[0]][_col[0]]) * cofactor00 +
      double(m[_row[0]][_col[1]]) * cofactor10 +
      double(m[_row[0]][_col[2]]) * cofactor20;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// ACCESS ALIGNED BOUNDING BOXES (AABB) AND BOUNDING VOLUME HIERARCHIES (BVH)
//
// These are used extensively in ray tracing.  See the hit() methods for 
// how tests are made against AABBs, BVH, and polygons.
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::AABB::AABB( const Model::real3& p )
{
    _min = p;
    _max = p;
}  

inline Model::AABB::AABB( const Model::real3& p0, const Model::real3& p1 )
{
    _min = p0;
    _max = p0;
    expand( p1 );
}  

inline Model::AABB::AABB( const Model::real3& p0, const Model::real3& p1, const Model::real3& p2 ) 
{
    _min = p0;
    _max = p0;
    expand( p1 );
    expand( p2 );
}  

inline void Model::AABB::pad( Model::real p ) 
{
    _min -= real3( p, p, p );
    _max += real3( p, p, p );
}

inline void Model::AABB::expand( const Model::AABB& other )
{
    for( uint i = 0; i < 3; i++ )
    {
        if ( other._min.c[i] < _min.c[i] ) _min.c[i] = other._min.c[i];
        if ( other._max.c[i] > _max.c[i] ) _max.c[i] = other._max.c[i];
    }
}

inline void Model::AABB::expand( const Model::real3& p ) 
{
    if ( p.c[0] < _min.c[0] ) _min.c[0] = p.c[0];
    if ( p.c[1] < _min.c[1] ) _min.c[1] = p.c[1];
    if ( p.c[2] < _min.c[2] ) _min.c[2] = p.c[2];
    if ( p.c[0] > _max.c[0] ) _max.c[0] = p.c[0];
    if ( p.c[1] > _max.c[1] ) _max.c[1] = p.c[1];
    if ( p.c[2] > _max.c[2] ) _max.c[2] = p.c[2];
}

inline bool Model::AABB::encloses( const real3& p ) const
{
    return _min.c[0] <= p.c[0] &&
           _min.c[1] <= p.c[1] &&
           _min.c[2] <= p.c[2] &&
           _max.c[0] >= p.c[0] &&
           _max.c[1] >= p.c[1] &&
           _max.c[2] >= p.c[2];
}

inline bool Model::AABB::encloses( const AABB& other ) const
{
    return _min.c[0] <= other._min.c[0] &&
           _min.c[1] <= other._min.c[1] &&
           _min.c[2] <= other._min.c[2] &&
           _max.c[0] >= other._max.c[0] &&
           _max.c[1] >= other._max.c[1] &&
           _max.c[2] >= other._max.c[2];
}

inline bool Model::AABB::overlaps( const Model::AABB& other ) const 
{
    return _min.c[0] <= other._max.c[0] &&
           _min.c[1] <= other._max.c[1] &&
           _min.c[2] <= other._max.c[2] &&
           _max.c[0] >= other._min.c[0] &&
           _max.c[1] >= other._min.c[1] &&
           _max.c[2] >= other._min.c[2];
}

static inline Model::real3 __transform_relative( const Model::real3& v, const Model::real3& center, const Model::real& sizeMaxRcp ) 
{ 
    return Model::real3( (v[0] - center[0]) * sizeMaxRcp + 0.5f, 
                         (v[1] - center[1]) * sizeMaxRcp + 0.5f, 
                         (v[2] - center[2]) * sizeMaxRcp + 0.5f );
}

inline Model::real3 Model::AABB::transform_relative( const Model::real3& v ) const
{ 
    return __transform_relative( v, center(), size_max_rcp() );
}

inline Model::AABB Model::AABB::transform_relative( const Model::AABB& other ) const
{
    real3 c = center(); 
    real  rcp = size_max_rcp();
    return AABB( __transform_relative( other._min, c, rcp ), __transform_relative( other._max, c, rcp ) );
}

inline bool Model::AABB::hit( const Model::real3& origin, const Model::real3& direction, const Model::real3& direction_inv, 
                              Model::real& t_min, Model::real& t_max ) const 
{
    mdout << "Model::AABB::hit: " << *this << " t_min=" << t_min << " t_max=" << t_max << "\n";
    (void)direction;
    for( uint a = 0; a < 3; a++ ) 
    {
        real dir_inv = direction_inv.c[a];
        real v0 = (_min.c[a] - origin.c[a]) * dir_inv;
        real v1 = (_max.c[a] - origin.c[a]) * dir_inv;
        t_min = std::fmax( t_min, std::fmin( v0, v1 ) );
        t_max = std::fmin( t_max, std::fmax( v0, v1 ) );
        mdout << "Model::AABB::hit:     " << a << ": _min=" << _min.c[a] << " _max=" << _max.c[a] << 
                                   " dir_inv=" << dir_inv << " origin=" << origin.c[a] << 
                                   " v0=" << v0 << " v1=" << v1 << " t_min=" << t_min << " t_max=" << t_max << "\n";
    }
    bool r = t_max >= std::fmax( t_min, real(0.0) );
    mdout << "Model::AABB::hit: return=" << r << "\n";
    return r;
}

inline Model::AABBD::AABBD( const Model::real3d& p )
{
    _min = p;
    _max = p;
}  

inline Model::AABBD::AABBD( const Model::real3d& p0, const Model::real3d& p1 )
{
    _min = p0;
    _max = p0;
    expand( p1 );
}  

inline Model::AABBD::AABBD( const Model::real3d& p0, const Model::real3d& p1, const Model::real3d& p2 ) 
{
    _min = p0;
    _max = p0;
    expand( p1 );
    expand( p2 );
}  

inline Model::AABBD::AABBD( const Model::AABB& box )
{
    _min = box._min.to_real3d();
    _max = box._max.to_real3d();
}  

inline void Model::AABBD::pad( Model::real64 p ) 
{
    _min -= real3d( p, p, p );
    _max += real3d( p, p, p );
}

inline void Model::AABBD::expand( const Model::AABBD& other )
{
    for( uint i = 0; i < 3; i++ )
    {
        if ( other._min.c[i] < _min.c[i] ) _min.c[i] = other._min.c[i];
        if ( other._max.c[i] > _max.c[i] ) _max.c[i] = other._max.c[i];
    }
}

inline void Model::AABBD::expand( const Model::real3d& p ) 
{
    if ( p.c[0] < _min.c[0] ) _min.c[0] = p.c[0];
    if ( p.c[1] < _min.c[1] ) _min.c[1] = p.c[1];
    if ( p.c[2] < _min.c[2] ) _min.c[2] = p.c[2];
    if ( p.c[0] > _max.c[0] ) _max.c[0] = p.c[0];
    if ( p.c[1] > _max.c[1] ) _max.c[1] = p.c[1];
    if ( p.c[2] > _max.c[2] ) _max.c[2] = p.c[2];
}

inline bool Model::AABBD::encloses( const real3d& p ) const
{
    return _min.c[0] <= p.c[0] &&
           _min.c[1] <= p.c[1] &&
           _min.c[2] <= p.c[2] &&
           _max.c[0] >= p.c[0] &&
           _max.c[1] >= p.c[1] &&
           _max.c[2] >= p.c[2];
}

inline bool Model::AABBD::encloses( const AABBD& other ) const
{
    return _min.c[0] <= other._min.c[0] &&
           _min.c[1] <= other._min.c[1] &&
           _min.c[2] <= other._min.c[2] &&
           _max.c[0] >= other._max.c[0] &&
           _max.c[1] >= other._max.c[1] &&
           _max.c[2] >= other._max.c[2];
}

inline bool Model::AABBD::overlaps( const Model::AABBD& other ) const 
{
    return _min.c[0] <= other._max.c[0] &&
           _min.c[1] <= other._max.c[1] &&
           _min.c[2] <= other._max.c[2] &&
           _max.c[0] >= other._min.c[0] &&
           _max.c[1] >= other._min.c[1] &&
           _max.c[2] >= other._min.c[2];
}

// Moeller code for fast triangle-box overlap test
// Taken from: https://github.com/3YOURMIND/3YDMoeller
//
inline bool axisTestX01(const Model::real3d &v0, const Model::real3d &v2, const Model::real3d &boxhalfsize,
                        const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                        Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p0 = a * v0[1] - b * v0[2];
    Model::real64 p2 = a * v2[1] - b * v2[2];
    if (p0 < p2)
    {
        min = p0;
        max = p2;
    }
    else
    {
        min = p2;
        max = p0;
    }
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool axisTestX2(const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &boxhalfsize,
                              const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                              Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p0 = a * v0[1] - b * v0[2];
    Model::real64 p1 = a * v1[1] - b * v1[2];
    if (p0 < p1)
    {
        min = p0;
        max = p1;
    }
    else
    {
        min = p1;
        max = p0;
    }
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool axisTestY02(const Model::real3d &v0, const Model::real3d &v2, const Model::real3d &boxhalfsize,
                               const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                               Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p0 = -a * v0[0] + b * v0[2];
    Model::real64 p2 = -a * v2[0] + b * v2[2];
    if (p0 < p2)
    {
        min = p0;
        max = p2;
    }
    else
    {
        min = p2;
        max = p0;
    }
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool axisTestY1(const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &boxhalfsize,
                              const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                              Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p0 = -a * v0[0] + b * v0[2];
    Model::real64 p1 = -a * v1[0] + b * v1[2];
    if (p0 < p1)
    {
        min = p0;
        max = p1;
    }
    else
    {
        min = p1;
        max = p0;
    }
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool axisTestZ12(const Model::real3d &v1, const Model::real3d &v2, const Model::real3d &boxhalfsize,
                               const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                               Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p1 = a * v1[0] - b * v1[1];
    Model::real64 p2 = a * v2[0] - b * v2[1];
    if (p2 < p1)
    {
        min = p2;
        max = p1;
    }
    else
    {
        min = p1;
        max = p2;
    }
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool axisTestZ0(const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &boxhalfsize,
                              const Model::real64 a, const Model::real64 b, const Model::real64 fa, const Model::real64 fb,
                              Model::real64 &min, Model::real64 &max, Model::real64 &rad)
{
    Model::real64 p0 = a * v0[0] - b * v0[1];
    Model::real64 p1 = a * v1[0] - b * v1[1];
    if (p0 < p1)
    {
        min = p0;
        max = p1;
    }
    else
    {
        min = p1;
        max = p0;
    }
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];
    if (min > rad || max < -rad)
    {
        return false;
    }
    return true;
}

inline bool planeBoxOverlap(const Model::real3d &normal, const Model::real3d &vert,
                            const Model::real3d &maxbox)  // -NJMP-
{
    size_t q;
    Model::real64 v;
    Model::real3d vmin, vmax;
    for (q = 0; q <= 2; q++)
    {
        v = vert[q];
        if (normal[q] > 0.0f)
        {
            vmin[q] = -maxbox[q] - v;
            vmax[q] = maxbox[q] - v;
        }
        else
        {
            vmin[q] = maxbox[q] - v;
            vmax[q] = -maxbox[q] - v;
        }
    }
    if (dot(normal, vmin) > 0.0f)
    {
        return false;
    }
    if (dot(normal, vmax) >= 0.0f)
    {
        return true;
    }
    return false;
}

bool triBoxOverlap( const Model::real3d &trivert0, const Model::real3d &trivert1, const Model::real3d &trivert2,
                    const Model::real3d &boxcenter, const Model::real3d &boxhalfsize )
{
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */

    Model::real3d v0, v1, v2;
    Model::real64 min, max, rad, fex, fey, fez;
    Model::real3d normal, e0, e1, e2;
    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    v0 = trivert0 - boxcenter;
    v1 = trivert1 - boxcenter;
    v2 = trivert2 - boxcenter;
    /* compute triangle edges */
    e0 = v1 - v0; /* tri edge 0 */
    e1 = v2 - v1; /* tri edge 1 */
    e2 = v0 - v2; /* tri edge 2 */

    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = std::fabs(e0[0]);
    fey = std::fabs(e0[1]);
    fez = std::fabs(e0[2]);

    if (!axisTestX01(v0, v2, boxhalfsize, e0[2], e0[1], fez, fey, min, max, rad))
    {
        return false;
    }
    if (!axisTestY02(v0, v2, boxhalfsize, e0[2], e0[0], fez, fex, min, max, rad))
    {
        return false;
    }
    if (!axisTestZ12(v1, v2, boxhalfsize, e0[1], e0[0], fey, fex, min, max, rad))
    {
        return false;
    }

    fex = std::fabs(e1[0]);
    fey = std::fabs(e1[1]);
    fez = std::fabs(e1[2]);

    if (!axisTestX01(v0, v2, boxhalfsize, e1[2], e1[1], fez, fey, min, max, rad))
    {
        return false;
    }
    if (!axisTestY02(v0, v2, boxhalfsize, e1[2], e1[0], fez, fex, min, max, rad))
    {
        return false;
    }
    if (!axisTestZ0(v0, v1, boxhalfsize, e1[1], e1[0], fey, fex, min, max, rad))
    {
        return false;
    }

    fex = std::fabs(e2[0]);
    fey = std::fabs(e2[1]);
    fez = std::fabs(e2[2]);

    if (!axisTestX2(v0, v1, boxhalfsize, e2[2], e2[1], fez, fey, min, max, rad))
    {
        return false;
    }
    if (!axisTestY1(v0, v1, boxhalfsize, e2[2], e2[0], fez, fex, min, max, rad))
    {
        return false;
    }
    if (!axisTestZ12(v1, v2, boxhalfsize, e2[1], e2[0], fey, fex, min, max, rad))
    {
        return false;
    }

    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */
    /* test in X-direction */

    minmax(v0[0], v1[0], v2[0], min, max);
    if (min > boxhalfsize[0] || max < -boxhalfsize[0])
    {
        return false;
    }
    /* test in Y-direction */
    minmax(v0[1], v1[1], v2[1], min, max);
    if (min > boxhalfsize[1] || max < -boxhalfsize[1])
    {
        return false;
    }
    /* test in Z-direction */
    minmax(v0[2], v1[2], v2[2], min, max);
    if (min > boxhalfsize[2] || max < -boxhalfsize[2])
    {
        return false;
    }

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    normal = cross(e0, e1);

    if (!planeBoxOverlap(normal, v0, boxhalfsize))
    {
        return false;
    }
    return true; /* box and triangle overlaps */
}

inline void reorderTriVertices( Model::real3d& p0, Model::real3d& p1, Model::real3d& p2, Model::real3d& w0, Model::real3d& w1, Model::real3d& w2, Model::real64& d0, Model::real64& d1, Model::real64& d2, uint32_t i )
{
    if ( i == 0 ) return;
    if ( i == 1 ) {
        Model::real3d pt = p0;
        Model::real3d wt = w0;
        Model::real64 dt = d0;
        p0 = p1; 
        w0 = w1;
        d0 = d1;
        p1 = p2;
        w1 = w2;
        d1 = d2;
        p2 = pt;
        w2 = wt;
        d2 = dt;
    } else {
        Model::real3d pt = p0;
        Model::real3d wt = w0;
        Model::real64 dt = d0;
        p0 = p2; 
        w0 = w2;
        d0 = d2;
        p2 = p1;
        w2 = w1;
        d2 = d1;
        p1 = pt;
        w1 = wt;
        d1 = dt;
    }
}

inline void applyTriWeights( const Model::real3d& v0, const Model::real3d& v1, const Model::real3d& v2, const Model::real3d& w, Model::real3d& p )
{
    Model::real64 w0 = w[0];
    Model::real64 w1 = w[1];
    Model::real64 w2 = 1.0 - w0 - w1;   // make sure they always add to 1.0
    assert( w0 >= 0.0 && w0 <= 1.0 );
    assert( w1 >= 0.0 && w1 <= 1.0 );
    assert( w2 >= 0.0 && w2 <= 1.0 );
    p = w0*v0 + w1*v1 + w2*v2;          // apply weights to original triangle coordinates
}

inline bool findTriPointInBox( const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &v2, 
                               const Model::AABBD& box, const Model::real3d& boxcenter, const Model::real3d& boxhalfsize, Model::real3d * p_ptr ) 
{
    //---------------------------------------------------------------------
    // This routine is called only if the triangle intersects the box, so
    // we can assume that. Make it so that vertex 0 is the closest to the center of the box.
    //---------------------------------------------------------------------
    Model::real3d p0 = v0;
    Model::real3d p1 = v1;
    Model::real3d p2 = v2;
    Model::real3d w0( 1.0, 0.0, 0.0 );  // w0[0]*v0 + w0[1]*v1 + w0[2]*v2
    Model::real3d w1( 0.0, 1.0, 0.0 );  
    Model::real3d w2( 0.0, 0.0, 1.0 );  
    Model::real3d tmp;
    mdout << "box=" << box << "\n";
    Model::real64 d0 = v0.distance_sqr( boxcenter ); 
    Model::real64 d1 = v1.distance_sqr( boxcenter ); 
    Model::real64 d2 = v2.distance_sqr( boxcenter ); 
    for( uint32_t i = 0;; i++ ) 
    {
        //---------------------------------------------------------------------
        // Reorder vertices by proximity to the center of the box.
        //---------------------------------------------------------------------
        reorderTriVertices( p0, p1, p2, w0, w1, w2, d0, d1, d2, (d0 <= d1 && d0 <= d2) ? 0 : (d1 <= d2) ? 1 : 2 );
        mdout << "    i=" << i << " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << " w0=" << w0 << " w1=" << w1 << " w2=" << w2 << " d0=" << d0 << " d1=" << d1 << " d2=" << d2 << "\n";

        //---------------------------------------------------------------------
        // See if any of the vertices are inside the box, favoring p0 which is the closest to the box center.
        //---------------------------------------------------------------------
        if ( box.encloses( p0 ) ) { *p_ptr = p0; mdout << "    chose p0\n"; return true; }
        if ( box.encloses( p1 ) ) { *p_ptr = p1; mdout << "    chose p1\n"; return true; }
        if ( box.encloses( p2 ) ) { *p_ptr = p2; mdout << "    chose p2\n"; return true; }

        //---------------------------------------------------------------------
        // Calculate midpoint of each side.
        // We do this to the weights, then apply the weights in a way that doesn't lose precision.
        //
        // We can use these midpoints to form 4 new triangles.
        // Switch to one that overlaps, favoring the p0-p01-p02 triangle, then middle triangle,
        // then one of the other two.
        //---------------------------------------------------------------------
        Model::real3d p01, p02, p12;
        Model::real3d w01 = (w0+w1).divby2(); 
        Model::real3d w02 = (w0+w2).divby2();
        Model::real3d w12 = (w1+w2).divby2(); 
        applyTriWeights( v0, v1, v2, w01, p01 );
        applyTriWeights( v0, v1, v2, w02, p02 );
        applyTriWeights( v0, v1, v2, w12, p12 );
        Model::real64 d01 = p01.distance_sqr( boxcenter );
        Model::real64 d02 = p02.distance_sqr( boxcenter );
        Model::real64 d12 = p12.distance_sqr( boxcenter );
        if ( triBoxOverlap( p0, p01, p02, boxcenter, boxhalfsize ) ) {
            mdout << "    switching to p0, p01, p02\n";
            p1 = p01;
            w1 = w01;
            d1 = d01;
            p2 = p02;
            w2 = w02;
            d2 = d02;
            continue;
        }
        if ( triBoxOverlap( p01, p12, p02, boxcenter, boxhalfsize ) ) {
            mdout << "    switching to p01, p12, p02\n";
            p0 = p01;
            w0 = w01;
            d0 = d01;
            p1 = p12;
            w1 = w12;
            d1 = d12;
            p2 = p02;
            w2 = w02;
            d2 = d02;
            continue;
        }
        if ( triBoxOverlap( p01, p1, p12, boxcenter, boxhalfsize ) ) {
            mdout << "    switching to p01, p1, p12\n";
            p0 = p01;
            w0 = w01;
            d0 = d01;
            p2 = p12;
            w2 = w12;
            d2 = d12;
            continue;
        }
        if ( triBoxOverlap( p02, p12, p2, boxcenter, boxhalfsize ) ) {
            mdout << "    switching to p02, p12, p2\n";
            p0 = p02;
            w0 = w02;
            d0 = d02;
            p1 = p12;
            w1 = w12;
            d1 = d12;
            continue;
        }
        std::cout << "ERROR: no subdivided triangle overlaps the box\n";
        assert( false );
    }
}

inline bool Model::AABBD::overlaps_triangle( const Model::real3d &v0, const Model::real3d &v1, const Model::real3d &v2, Model::real3d * p_ptr ) const
{
    real3d boxcenter   = (_min + _max) * 0.5;  
    real3d boxhalfsize = (_max - _min) * 0.5;
    if ( !triBoxOverlap( v0, v1, v2, boxcenter, boxhalfsize ) ) return false;
    if ( p_ptr != nullptr ) {
        //---------------------------------------------------------------------
        // For now, use recursive triangle subdivision to find any (sub)triangle center point
        // that is inside this box.
        //---------------------------------------------------------------------
        assert( findTriPointInBox( v0, v1, v2, *this, boxcenter, boxhalfsize, p_ptr ) ); 
    }
    return true;
}

static inline Model::real3d __transform_relative( const Model::real3d& v, const Model::real3d& center, const Model::real64& sizeMaxRcp ) 
{ 
    return Model::real3d( (v[0] - center[0]) * sizeMaxRcp + 0.5f, 
                          (v[1] - center[1]) * sizeMaxRcp + 0.5f, 
                          (v[2] - center[2]) * sizeMaxRcp + 0.5f );
}

inline Model::real3d Model::AABBD::transform_relative( const Model::real3d& v ) const
{ 
    return __transform_relative( v, center(), size_max_rcp() );
}

inline Model::AABBD Model::AABBD::transform_relative( const Model::AABBD& other ) const
{
    real3d c = center(); 
    real64 rcp = size_max_rcp();
    return AABBD( __transform_relative( other._min, c, rcp ), __transform_relative( other._max, c, rcp ) );
}

inline bool Model::AABBD::hit( const Model::real3d& origin, const Model::real3d& direction, const Model::real3d& direction_inv, 
                               Model::real64& t_min, Model::real64& t_max ) const 
{
    mdout << "Model::AABBD::hit: " << *this << " t_min=" << t_min << " t_max=" << t_max << "\n";
    (void)direction;
    for( uint a = 0; a < 3; a++ ) 
    {
        real64 dir_inv = direction_inv.c[a];
        real64 v0 = (_min.c[a] - origin.c[a]) * dir_inv;
        real64 v1 = (_max.c[a] - origin.c[a]) * dir_inv;
        t_min = std::fmax( t_min, std::fmin( v0, v1 ) );
        t_max = std::fmin( t_max, std::fmax( v0, v1 ) );
        mdout << "Model::AABBD::hit:     " << a << ": _min=" << _min.c[a] << " _max=" << _max.c[a] << 
                                   " dir_inv=" << dir_inv << " origin=" << origin.c[a] << 
                                   " v0=" << v0 << " v1=" << v1 << " t_min=" << t_min << " t_max=" << t_max << "\n";
    }
    bool r = t_max >= std::fmax( t_min, real64(0.0) );
    mdout << "Model::AABBD::hit: return=" << r << "\n";
    return r;
}

inline Model::AABBI::AABBI( const Model::_int p[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p[i];
        _max[i] = p[i];
    }
}  

inline Model::AABBI::AABBI( const Model::_int p0[], const _int p1[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p0[i];
        _max[i] = p0[i];
    }
    expand( p1 );
}  

inline Model::AABBI::AABBI( const Model::_int p0[], const _int p1[], const _int p2[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p0[i];
        _max[i] = p0[i];
    }
    expand( p1 );
    expand( p2 );
}  

inline void Model::AABBI::expand( const Model::_int p[] ) 
{
    if ( p[0] < _min[0] ) _min[0] = p[0];
    if ( p[1] < _min[1] ) _min[1] = p[1];
    if ( p[2] < _min[2] ) _min[2] = p[2];
    if ( p[0] > _max[0] ) _max[0] = p[0];
    if ( p[1] > _max[1] ) _max[1] = p[1];
    if ( p[2] > _max[2] ) _max[2] = p[2];
}

inline bool Model::AABBI::encloses( const Model::_int p[] ) const
{
    return _min[0] <= p[0] &&
           _min[1] <= p[1] &&
           _min[2] <= p[2] &&
           _max[0] >= p[0] &&
           _max[1] >= p[1] && 
           _max[2] >= p[2];
}

inline bool Model::AABBI::encloses( Model::_int x, Model::_int y, Model::_int z ) const
{
    const _int p[3] = { x, y, z };
    return encloses( p );
}

inline Model::AABBU64::AABBU64( const Model::uint64 p[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p[i];
        _max[i] = p[i];
    }
}  

inline Model::AABBU64::AABBU64( const Model::uint64 p0[], const uint64 p1[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p0[i];
        _max[i] = p0[i];
    }
    expand( p1 );
}  

inline Model::AABBU64::AABBU64( const Model::uint64 p0[], const uint64 p1[], const uint64 p2[] )
{
    for( uint32_t i = 0; i < 3; i++ ) 
    {
        _min[i] = p0[i];
        _max[i] = p0[i];
    }
    expand( p1 );
    expand( p2 );
}  

inline void Model::AABBU64::expand( const Model::uint64 p[] ) 
{
    if ( p[0] < _min[0] ) _min[0] = p[0];
    if ( p[1] < _min[1] ) _min[1] = p[1];
    if ( p[2] < _min[2] ) _min[2] = p[2];
    if ( p[0] > _max[0] ) _max[0] = p[0];
    if ( p[1] > _max[1] ) _max[1] = p[1];
    if ( p[2] > _max[2] ) _max[2] = p[2];
}

inline bool Model::AABBU64::encloses( const Model::uint64 p[] ) const
{
    return _min[0] <= p[0] &&
           _min[1] <= p[1] &&
           _min[2] <= p[2] &&
           _max[0] >= p[0] &&
           _max[1] >= p[1] && 
           _max[2] >= p[2];
}

inline bool Model::AABBU64::encloses( Model::uint64 x, Model::uint64 y, Model::uint64 z ) const
{
    const uint64 p[3] = { x, y, z };
    return encloses( p );
}

inline Model::Ray::Ray( const Model::real3& a, const Model::real3& b, Model::RAY_KIND kind, Model::real32 solid_angle, Model::real32 cone_angle ) 
        : A(a), B(b), KIND(kind), SOLID_ANGLE(solid_angle), CONE_ANGLE(cone_angle) 
{
    B_norm = B; 
    B_norm.normalize(); 
    B_inv = Model::real3(1,1,1); 
    B_inv /= B; 
    top = -1;
}

inline void Model::Ray::init_normalized( const Model::real3& a, const Model::real3& b, RAY_KIND kind, Model::real32 solid_angle, Model::real32 cone_angle ) 
{
    A = a; 
    B = b; 
    SOLID_ANGLE = solid_angle;
    CONE_ANGLE = cone_angle;
    B_norm = b; 
    B_inv = Model::real3(1,1,1); 
    B_inv /= B; 
    KIND = kind;
    top = -1;
}

inline Model::real3 Model::Ray::point_at_parameter( Model::real t ) const 
{ 
    return A + t*B; 
}

inline Model::real3 Model::Ray::normalized_point_at_parameter( Model::real t ) const 
{ 
    return A + t*B_norm; 
}

inline const Model::real3 Model::Ray::cone_base_dir_x( void ) const 
{ 
    // Get the cone edge vector in 2D by rotating direction by CONE_ANGLE.
    //
    real32 sin, cos;
    sincos(CONE_ANGLE, sin, cos);
    real32 x = direction()[0];
    real32 y = direction()[1];
    real32 rx = x*cos - y*sin;
    real32 ry = x*sin + y*cos;

    // base_dir_x = edge - direction
    //
    return real3(rx-x, ry-y, 0);
}

inline Model::Ray64::Ray64( const Model::real3d& a, const Model::real3d& b, Model::RAY_KIND kind, Model::real64 solid_angle, Model::real64 cone_angle ) 
        : A(a), B(b), KIND(kind), SOLID_ANGLE(solid_angle), CONE_ANGLE(cone_angle) 
{
    B_norm = B; 
    B_norm.normalize(); 
    B_inv = Model::real3d(1,1,1); 
    B_inv /= B; 
    top = -1;
}

inline void Model::Ray64::init_normalized( const Model::real3d& a, const Model::real3d& b, RAY_KIND kind, Model::real64 solid_angle, Model::real64 cone_angle ) 
{
    A = a; 
    B = b; 
    SOLID_ANGLE = solid_angle;
    CONE_ANGLE = cone_angle;
    B_norm = b; 
    B_inv = Model::real3d(1,1,1); 
    B_inv /= B; 
    KIND = kind;
    top = -1;
}

inline Model::real3d Model::Ray64::point_at_parameter( Model::real64 t ) const 
{ 
    return A + t*B; 
}

inline Model::real3d Model::Ray64::normalized_point_at_parameter( Model::real64 t ) const 
{ 
    return A + t*B_norm; 
}

inline const Model::real3d Model::Ray64::cone_base_dir_x( void ) const 
{ 
    // Get the cone edge vector in 2D by rotating direction by CONE_ANGLE.
    //
    real64 sin, cos;
    sincos(CONE_ANGLE, sin, cos);
    real64 x = direction()[0];
    real64 y = direction()[1];
    real64 rx = x*cos - y*sin;
    real64 ry = x*sin + y*cos;

    // base_dir_x = edge - direction
    //
    return real3d(rx-x, ry-y, 0);
}

Model::Camera::Camera( real3 lookfrom, real3 lookat, real3 vup, real vfov, real near, real far, real aspect, real aperture, real focus_dist, int nx, int ny, int spp, uint name_i )
    : lookfrom(lookfrom), lookat(lookat), vup(vup), aperture(aperture), focus_dist(focus_dist), near(near), far(far), vfov(vfov), aspect(aspect), name_i(name_i)
{ 
    lens_radius = divby2(aperture);
    real theta = vfov*PI_DIV_180;
    real half_height = std::tan(divby2(theta));
    real half_width = aspect * half_height;
    w = (lookfrom - lookat).normalized();
    u = cross(vup, w).normalized();
    v = cross(w, u);
    lower_left_corner = lookfrom - half_width*focus_dist*u -half_height*focus_dist*v - focus_dist*w;
    horizontal = mulby2(half_width*focus_dist)*u;
    vertical = mulby2(half_height*focus_dist)*v;
    real32 area = cross(horizontal, vertical).length();
    real3 screen_center = lower_left_corner + (horizontal + vertical).divby2();
    real32 distance_to_screen = (lookfrom - screen_center).length();
    solid_angle = area / (distance_to_screen*real32(nx*ny*spp));

    //------------------------------------------------
    // Construct the perspective projection matrix.
    // perspective_projection = perspective * view
    //------------------------------------------------
    Matrix V = Matrix::make_view( lookfrom, lookat, vup );
    Matrix P = Matrix::make_perspective( vfov, aspect, near, far );
    P.transform( V, perspective_projection );

    //------------------------------------------------
    // Extract frustum_planes from perspective_projection matrix.
    // https://www.gamedevs.org/uploads/fast-extraction-viewing-frustum-planes-from-world-view-projection-matrix.pdf
    // 0=near 1=far 2=bottom 3=top left=4 right=5
    //------------------------------------------------
    const Matrix& pp = perspective_projection;
    frustum_planes[0].normal   = real3( pp.m[3][0] + pp.m[2][0], pp.m[3][1] + pp.m[2][1], pp.m[3][2] + pp.m[2][2] );
    frustum_planes[0].distance = pp.m[3][3] + pp.m[2][3];
    frustum_planes[1].normal   = real3( pp.m[3][0] - pp.m[2][0], pp.m[3][1] - pp.m[2][1], pp.m[3][2] - pp.m[2][2] );
    frustum_planes[1].distance = pp.m[3][3] - pp.m[2][3];
    frustum_planes[2].normal   = real3( pp.m[3][0] + pp.m[1][0], pp.m[3][1] + pp.m[1][1], pp.m[3][2] + pp.m[1][2] );
    frustum_planes[2].distance = pp.m[3][3] + pp.m[1][3];
    frustum_planes[3].normal   = real3( pp.m[3][0] - pp.m[1][0], pp.m[3][1] - pp.m[1][1], pp.m[3][2] - pp.m[1][2] );
    frustum_planes[3].distance = pp.m[3][3] - pp.m[1][3];
    frustum_planes[4].normal   = real3( pp.m[3][0] + pp.m[0][0], pp.m[3][1] + pp.m[0][1], pp.m[3][2] + pp.m[0][2] );
    frustum_planes[4].distance = pp.m[3][3] + pp.m[0][3];
    frustum_planes[5].normal   = real3( pp.m[3][0] - pp.m[0][0], pp.m[3][1] - pp.m[0][1], pp.m[3][2] - pp.m[0][2] );
    frustum_planes[5].distance = pp.m[3][3] - pp.m[0][3];
}

inline void Model::Camera::concentric_point_on_unit_disk( real s, real t, real& x, real& y ) const
{
    // more complex than polar, but according to 
    // A Realistic Camera Model for Computer Graphics (Kolb et al)
    // has 15% improvement in error over polar for camera lens sampling
    // possibly similar improvement for sampling cones?
    real a = mulby2(s) - 1.0;
    real b = mulby2(t) - 1.0;
    real rad, theta;

    if (std::abs(a) > std::abs(b)) {
        rad = a;
        theta = PI_DIV_4*(b/a);
    } else {
        rad = b;
        theta = PI_DIV_2 - PI_DIV_4*(a/b);
    }
    sincos(theta, y, x, rad);
}

Model::Ray Model::Camera::get_ray( real s_pixel, real t_pixel, real s_lens, real t_lens ) const 
{
    real _x, _y;
    concentric_point_on_unit_disk(s_lens, t_lens, _x, _y);
    real32 x = _x;
    real32 y = _y;
    x *= lens_radius;
    y *= lens_radius;
    real3 offset;
    for(int i = 0; i < 3; i++)
    {
        offset[i] = real32(u[0]) * x + real32(v[0]) * y;
    }
    Ray r = Ray( lookfrom + offset, lower_left_corner + s_pixel*horizontal + t_pixel*vertical - lookfrom - offset, RAY_KIND::DIRECT, solid_angle ); 
    mdout << "camera::get_ray: s_lens=" << s_lens << " t_lens=" << t_lens << " x=" << x << " y=" << y << 
                            " u=" << u << " v=" << v << " lens_radius=" << lens_radius << " offset=" << offset <<
                            " lower_left_corner=" << lower_left_corner << " s_pixel=" << s_pixel << " t_pixel=" << t_pixel <<
                            " solid_angle=" << solid_angle << " r=" << r << "\n";
    return r;
}

// see: "Optimized View Frustum Culling Algorithms", http://www.ce.chalmers.se/staff/uffe
//
bool Model::Camera::overlaps_box( const Model::AABB& box ) const 
{
    real3 n, p;
    for( uint32_t i = 0; i < 6; i++ )
    {
        p = box._min;
        n = box._max;

        for( uint32_t j = 0; j < 3; j++ ) 
        {
            if ( frustum_planes[i].normal[j] >= 0.0 ) {
                p[j] = box._max[j];
                n[j] = box._min[j];
            }
        }

        real dot = n.dot( frustum_planes[i].normal ) + frustum_planes[i].distance;
        if ( dot > 0.0 ) return false;
    }

    return true;
}

bool Model::Polygon::bounding_box( const Model * model, Model::AABB& box, real padding ) const 
{
    const Vertex * vertexes = &model->vertexes[vtx_i];
    for( uint32_t i = 0; i < vtx_cnt; i++ )
    {
        const real3& p = model->positions[vertexes[i].v_i];
        if ( i == 0 ) {
            box = AABB( p );
        } else {
            box.expand( p );
        }
    }
    box.pad( padding );
    return true;
}

bool Model::Polygon::vertex_info( const Model * model, uint _vtx_cnt, real3 p[], real3 n[], real2 uv[] ) const
{
    if (_vtx_cnt != vtx_cnt) return false;
    for( uint i = 0; i < vtx_cnt; i++ ) 
    {
        const Vertex * vertex = &model->vertexes[vtx_i+i];
        p[i] = model->positions[vertex->v_i];
        if ( vertex->vn_i != uint(-1) ) {
            n[i] = model->normals[vertex->vn_i];
        } else {
            n[i] = normal;
        }
        if ( vertex->vt_i != uint(-1) ) {
            uv[i] = model->texcoords[vertex->vt_i];
        } else {
            uv[i] = real2(0,0);
        }
    }
    return true;
}

bool Model::Polygon::vertex_info( const Model * model, uint _vtx_cnt, real3 p[] ) const
{
    if (_vtx_cnt != vtx_cnt) return false;
    for( uint i = 0; i < vtx_cnt; i++ ) 
    {
        const Vertex * vertex = &model->vertexes[vtx_i+i];
        p[i] = model->positions[vertex->v_i];
    }
    return true;
}

bool Model::Polygon::sample( const Model * model, HitRecord& rec, real solid_angle, real dir_sqr ) const
{
    if ( vtx_cnt != 3 ) return false;

    uint poly_i = this - model->polygons;
    const Vertex * vertexes = &model->vertexes[vtx_i];
    const real3 * positions = model->positions;
    const real3& p0 = positions[vertexes[0].v_i];
    const real3& p1 = positions[vertexes[1].v_i];
    const real3& p2 = positions[vertexes[2].v_i];
    const real3 p_m_p0  = rec.p - p0;
    const real3 p1_m_p0 = p1    - p0;
    const real3 p2_m_p0 = p2    - p0;
    real area1 = p1_m_p0.cross( p_m_p0 ).dot( normal );
    real area2 = p_m_p0.cross( p2_m_p0 ).dot( normal );
    real area_2x = area * 2.0;
    real beta    = area1/area_2x;
    real gamma   = area2/area_2x;
    mdout << "Model::Polygon::sample: poly_i=" << poly_i << " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << " beta=" << beta << " gamma=" << gamma << "\n";
    if ( beta < 0.0 || gamma < 0.0 || (beta + gamma) > 1.0 ) {
        mdout << "Model::Polygon::sample: poly_i=" << poly_i << " beta=" << beta << " gamma=" << gamma << ", so REJECTED\n";
        return false;
    }
    real alpha = 1.0 - beta - gamma;

    const real3 * normals   = model->normals;
    const real2 * texcoords = model->texcoords;
    const real3& n0 = normals[vertexes[0].vn_i];
    const real3& n1 = normals[vertexes[1].vn_i];
    const real3& n2 = normals[vertexes[2].vn_i];
    real u0, u1, u2, v0, v1, v2;
    if ( vertexes[0].vt_i != uint(-1) ) {
        u0 = texcoords[vertexes[0].vt_i].c[0];
        u1 = texcoords[vertexes[1].vt_i].c[0];
        u2 = texcoords[vertexes[2].vt_i].c[0];
        v0 = texcoords[vertexes[0].vt_i].c[1];
        v1 = texcoords[vertexes[1].vt_i].c[1];
        v2 = texcoords[vertexes[2].vt_i].c[1];
    } else {
        u0 = 0.0;
        u1 = 0.0;
        u2 = 0.0;
        v0 = 0.0;
        v1 = 0.0;
        v2 = 0.0;
    }

    if ( solid_angle != 0.0 ) {
        real distance_squared = rec.t*rec.t / dir_sqr;
        real ray_footprint_area_on_triangle = solid_angle*distance_squared;
        real uv_area_of_triangle = divby4(std::abs(u0*v1 + u1*v2 + u2*v0 - u0*v2 - u1*v0 - u2*v1));
        rec.frac_uv_cov = ray_footprint_area_on_triangle * uv_area_of_triangle / area;
    } else {
        rec.frac_uv_cov = 0.0;
    }

    rec.poly_i = poly_i;
    rec.grid_i = uint(-1);
    rec.normal = (std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2])) ? real3(0,1,0) : normal;
    if ( vertexes[0].vn_i != uint(-1) ) {
        rec.shading_normal = n0*alpha + n1*gamma + n2*beta;
        rec.shading_normal.normalize();
        if (std::isnan(rec.shading_normal[0]) || std::isnan(rec.shading_normal[1]) || std::isnan(rec.shading_normal[2])) rec.shading_normal = rec.normal;
    } else {
        rec.shading_normal = normal;
    }

    rec.u = alpha*u0 + gamma*u1 + beta*u2 ;
    rec.v = alpha*v0 + gamma*v1 + beta*v2 ;

    real3 deltaPos1 = p1-p0;
    real3 deltaPos2 = p2-p0;

    real deltaU1 = u1-u0;
    real deltaU2 = u2-u0;
    real deltaV1 = v1-v0;
    real deltaV2 = v2-v0;

    real r = 1.0 / (deltaU1 * deltaV2 - deltaV1 * deltaU2);
    rec.tangent = (deltaPos1 * deltaV2 - deltaPos2 * deltaV1)*r;
    rec.tangent_normalized = rec.tangent.normalize();
    rec.bitangent = (deltaPos2 * deltaU1 - deltaPos1 * deltaU2)*r;

    rec.model = model;
    return true;
}

bool Model::Polygon::hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv,
        real solid_angle, real t_min, real t_max, HitRecord& rec ) const
{
    (void)direction_inv;
    if ( vtx_cnt != 3 ) return false;

    // triangle - this code comes originally from Peter Shirley
    const Vertex * vertexes = &model->vertexes[vtx_i];
    const real3 * positions = model->positions;
    const real3& p0 = positions[vertexes[0].v_i];

    // plane equation (p - corner) dot N = 0
    // (o + t*v - corner) dot N = 0
    // t*dot(v,N) = (corner-o) dot N
    // t = dot(corner - o, N) / dot(v,N)
    real d = direction.dot( normal );
    real t = (p0 - origin).dot( normal ) / d;
    uint poly_i = this - model->polygons;
    mdout << "Model::Polygon::hit: model=" << model << " poly_i=" << poly_i << " origin=" << origin << 
                                 " direction=" << direction << " direction_inv=" << direction_inv << " solid_angle=" << solid_angle <<
                                 " d=" << d << " t=" << t << " t_min=" << t_min << " t_max=" << t_max << "\n";
    if ( t <= t_min || t >= t_max ) {
        mdout << "Model::Polygon::hit: NOT a hit, poly_i=" << poly_i << "\n";
        return false;
    }

    rec.t = t;
    rec.p = origin + direction * t;
    if ( !sample( model, rec, solid_angle, direction.length_sqr() ) ) {
        mdout << "Model::Polygon::hit: NOT a hit, poly_i=" << poly_i << "\n";
        return false;
    }

    if (mtl_i != uint(-1) && model->materials[mtl_i].map_d_i != uint(-1)) {
        // code slightly modified from texture.h   We only need one compoent
        // From there, the first texel will be at &texels[texel_i]. 
        // Pull out the map_d_i.  If it's -1, there's no alpha texture.  Also use that to get to the Texture using model->textures[map_d_i].
        uint texel_i = model->textures[model->materials[mtl_i].map_d_i].texel_i; 
        unsigned char *mdata = &(model->texels[texel_i]);
        int nx = model->textures[model->materials[mtl_i].map_d_i].width;
        int ny = model->textures[model->materials[mtl_i].map_d_i].height;
        int mx = nx;
        int my = ny;
        real u = rec.u;
        real v = rec.v;
        real sqrt_nx_ny = std::sqrt(nx*ny);
        real width_of_footprint = std::sqrt(rec.frac_uv_cov) * sqrt_nx_ny;
        real mip_level = std::log2( width_of_footprint );
        int nchan = model->textures->nchan;
        for (int imip_level = mip_level; imip_level > 0 && !(mx <= 1 && my <= 1); imip_level--)
        {
            // find the proper mip texture
            mdata += nchan * mx * my;
            if ( mx != 1 ) mx >>= 1;
            if ( my != 1 ) my >>= 1;
        }
        if (std::isnan(u)) {
            u = 0.0;
        }
        if (std::isnan(v)) {
            v = 0.0;
        }
        if (u < 0.0) {
            int64_t i = u;
            u -= i-1;
        }
        if (v < 0.0) {
            int64_t i = v;
            v -= i-1;
        }
        if (u >= 1.0) {
            int64_t i = u;
            u -= i;
        }
        if (v >= 1.0) {
            int64_t i = v;
            v -= i;
        }
        int i = (    u)*real(mx);
        int j = (1.0-v)*real(my);
        while (i >= mx && mx != 0) i -= mx;
        while (j >= my && my != 0) j -= my;

        float opacity = float(mdata[nchan*i + nchan*mx*j+3]) / 255.0;

        if (float(MODEL_UNIFORM_FN()) > opacity) {
            mdout << "Model::Polygon::hit: NOT HIT poly_i=" << poly_i << " uniform() > opacity=" << opacity << "\n";
            return false;
        }
    }

    mdout << "Model::Polygon::hit: HIT poly_i=" << poly_i << " t=" << rec.t << 
             " p=" << rec.p << " normal=" << rec.normal << 
             " frac_uv_cov=" << rec.frac_uv_cov << 
             " u=" << rec.u << " v=" << rec.v << " mtl_i=" << mtl_i << "\n";
    return true;
}

// simple standalone triangle hit that implements a small subset of the previous for generic usage
//
inline bool triangle_hit( const Model::real3& origin, const Model::real3& direction, Model::real t_min, Model::real t_max,
                          const Model::real3& p0, const Model::real3& p1, const Model::real3& p2, Model::real& t, Model::real3& p )
{
    Model::real3 normal = (p1 - p0).cross( p2 - p0 );
    Model::real d = direction.dot( normal );
    t = (p0 - origin).dot( normal ) / d;
    if ( t <= t_min || t >= t_max ) {
        return false;
    } else {
        p = origin + direction * t;
        return true;
    }
}

inline bool triangle_hit( const Model::real3d& origin, const Model::real3d& direction, Model::real64 t_min, Model::real64 t_max,
                          const Model::real3d& p0, const Model::real3d& p1, const Model::real3d& p2, Model::real64& t, Model::real3d& p )
{
    Model::real3d normal = (p1 - p0).cross( p2 - p0 );
    Model::real64 d = direction.dot( normal );
    t = (p0 - origin).dot( normal ) / d;
    if ( t <= t_min || t >= t_max ) {
        return false;
    } else {
        p = origin + direction * t;
        return true;
    }
}

bool Model::Polygon::is_emissive( const Model * model ) const
{
    return mtl_i != uint(-1) && model->materials[mtl_i].is_emissive();
}

bool Model::Volume::bounding_box( const Model * model, AABB& box, real padding ) const
{
    // all grids in a volume must have the same bbox
    //
    return model->volume_grids[0].bounding_box( model, box, padding );
}

bool Model::VolumeGrid::bounding_box( const Model * model, Model::AABB& box, real padding ) const 
{
    (void)model;
    box = world_box;
    box.pad( padding );
    return true;
}

bool Model::VolumeGrid::bounding_box( const Model * model, _int x, _int y, _int z, AABBD& box ) const
{
    (void)model;
    real3d xyz_r( x, y, z );
    box._min = xyz_r * world_voxel_size;
    box._max = box._min + world_voxel_size;
    return true;
}

const void * Model::VolumeGrid::value_ptr( const Model * model, _int x, _int y, _int z, bool * is_active ) const 
{
    //-------------------------------------------------------------------
    // This is where the crazy stuff happens.  
    // We walk the tree and attempt to find the voxel value at [x, y, z].
    // If the value is inactive, then we return nullptr.
    //
    // All the tree data structures are localized to this routine.
    //
    // Some of the data structures are variable-length, so we split them
    // before and after variable-length fields.
    //-------------------------------------------------------------------
    mdout << "VOLUME: xyz=[" << x << "," << y << "," << z << "]\n";

    // affine transform and its inverse represented as a 3x3 matrix and a real3 translation
    //
    struct MapData { 
        float  mMatF[9]; 
        float  mInvMatF[9];
        float  mVecF[3];
        float  mTaperF;
        double mMatD[9];
        double mInvMatD[9];
        double mVecD[3];
        double mTaperD;
    };

    // constants for 4-level tree structure
    //
    constexpr uint LEAF_LEVEL         = 0;
    constexpr uint LEAF_LOG2DIM       = 3;
    constexpr uint LEAF_TOTAL         = LEAF_LOG2DIM; // dimensions in index space
    constexpr uint LEAF_SIZE          = 1 << (3 * LEAF_LOG2DIM);
    constexpr uint LEAF_MASK_SIZE     = LEAF_SIZE >> 3;                       
    constexpr uint LEAF_MASK          = (1 << LEAF_TOTAL) - 1; 

    constexpr uint INTERNAL1_LEVEL    = LEAF_LEVEL + 1;
    constexpr uint INTERNAL1_LOG2DIM  = LEAF_LOG2DIM + 1;
    constexpr uint INTERNAL1_TOTAL    = INTERNAL1_LOG2DIM + LEAF_TOTAL;            // dimensions in index space
    constexpr uint INTERNAL1_SIZE     = 1 << (3 * INTERNAL1_LOG2DIM);              // number of tile values
    constexpr uint INTERNAL1_MASK_SIZE= INTERNAL1_SIZE >> 3;                       // bytes in value or child mask
    constexpr uint INTERNAL1_MASK     = (1 << INTERNAL1_TOTAL) - 1; 

    constexpr uint INTERNAL2_LEVEL    = INTERNAL1_LEVEL + 1;
    constexpr uint INTERNAL2_LOG2DIM  = INTERNAL1_LOG2DIM + 1;
    constexpr uint INTERNAL2_TOTAL    = INTERNAL2_LOG2DIM + INTERNAL1_TOTAL;
    constexpr uint INTERNAL2_SIZE     = 1 << (3 * INTERNAL2_LOG2DIM);
    constexpr uint INTERNAL2_MASK_SIZE= INTERNAL2_SIZE >> 3;                       
    constexpr uint INTERNAL2_MASK     = (1 << INTERNAL2_TOTAL) - 1; 

    constexpr uint ROOT_LEVEL         = INTERNAL2_LEVEL + 1;      

    //-------------------------------------------------------------------
    // The value size doesn't change.
    // Figure that out now.
    //-------------------------------------------------------------------
    uint value_size;
    switch( voxel_type ) 
    {
        case VolumeVoxelType::FLOAT:    value_size = sizeof(float);     break;
        case VolumeVoxelType::DOUBLE:   value_size = sizeof(double);    break;
        case VolumeVoxelType::INT16:    value_size = sizeof(int16_t);   break;
        case VolumeVoxelType::INT32:    value_size = sizeof(int32_t);   break;
        case VolumeVoxelType::INT64:    value_size = sizeof(int64_t);   break;
        case VolumeVoxelType::VEC3F:    value_size = 3*sizeof(float);   break;
        case VolumeVoxelType::VEC3D:    value_size = 3*sizeof(double);  break;
        case VolumeVoxelType::MASK:     value_size = sizeof(uint64_t);  break;
        case VolumeVoxelType::FP16:     value_size = 2;                 break;
        default:                        value_size = 0; die_assert( false, "bad voxel_type" ); break;
    }

    //-------------------------------------------------------------------
    // The GridData lives starts at &voxels[voxel_i].
    // TreeData comes after the name.
    // RootData comes after the TreeData.
    //-------------------------------------------------------------------
    static const int MaxNameSize = 256;
    struct GridData { 
        uint64_t    mMagic;             // 8 byte magic to validate it is valid grid data.
                                        // REFERENCEMAGIC = 0x4244566f6e614eL;
        char        mGridName[MaxNameSize];
        real64      mBBoxMin[3];        // min corner of AABB
        real64      mBBoxMax[3];        // max corner of AABB
        MapData     mMap;               // affine transformation between index and world space in both single and double precision
        double      mUniformScale;      // size of a voxel in world units
        VolumeGridClass mGridClass;     // 4 bytes
        VolumeVoxelType mVoxelType;     // 4 bytes
        uint32_t    mBlindDataCount;    // count of GridBlindMetaData structures that follow this grid (after the gridname).
        uint32_t    padding[2];         // seems to be needed to match up with NanoVDB
    };

    const uint8_t * grid_data = &model->voxels[voxel_i];
    const GridData * grid = reinterpret_cast< const GridData * >( grid_data );
    die_assert( grid->mMagic == 0x304244566f6e614eUL, "bad magic number in GridData" );
    die_assert( grid->mBlindDataCount == 0, "mBlindDataCount must be 0 currently" );
    mdout << "VOLUME: grid=" << grid << " mUniformScale=" << "\n";

    struct TreeData
    {
        uint64_t mBytes[ROOT_LEVEL + 1]; // byte offsets to nodes of each type
        uint32_t mCount[ROOT_LEVEL + 1]; // total number of nodes of each type
    };

    const uint8_t * tree_data = grid_data + sizeof(GridData);
    const TreeData * tree = reinterpret_cast< const TreeData * >( tree_data );
    mdout << "VOLUME: tree=" << tree << " mBytes[ROOT_LEVEL]=" << tree->mBytes[ROOT_LEVEL] << "\n";

    struct RootData4
    {
        float    mBBoxMin[3];           // min corner of world bbox
        float    mBBoxMax[3];           // max corner of world bbox
        uint64_t mActiveVoxelCount;     // total number of active voxels in the root and all its child nodes
        uint32_t mTileCount;            // number of tiles and child pointers in the root node
        uint32_t mPadding1[3];          // to match reference
        float    mBackGround;           // background value, i.e., value of any unset voxel
        float    mValueMin;             // minimum value
        float    mValueMax;             // maximum value
        uint32_t mPadding2[1];          // to match reference
    };

    struct RootData12
    {
        float    mBBoxMin[3];           // min corner of world bbox
        float    mBBoxMax[3];           // max corner of world bbox
        uint64_t mActiveVoxelCount;     // total number of active voxels in the root and all its child nodes
        uint32_t mTileCount;            // number of tiles and child pointers in the root node
        uint32_t mPadding1[3];          // to match reference
        float    mBackGround[3];        // background value, i.e., value of any unset voxel
        float    mValueMin[3];          // minimum value
        float    mValueMax[3];          // maximum value
        uint32_t mPadding2[3];          // to match reference
    };

    die_assert( value_size == 4 || value_size == 12, "unexpected value_size" );

    const uint8_t *   root_data = tree_data + tree->mBytes[ROOT_LEVEL];
    const uint        root_size = (value_size == 4) ? sizeof(RootData4) : sizeof(RootData12);
    const RootData4 * root4  = reinterpret_cast< const RootData4 *>( root_data );
    const uint8_t *   root_tile_data = reinterpret_cast< const uint8_t * >( root_data + root_size );
    mdout << "VOLUME: root=" << reinterpret_cast<const void *>(root_data) << " root_size=" << root_size << 
            " &mTileCount=" << reinterpret_cast<const void *>(&root4->mTileCount) << 
            " &mActiveVoxelCount=" << reinterpret_cast<const void *>(&root4->mActiveVoxelCount) << 
            " mTileCount=" << root4->mTileCount << " mActiveVoxelCount=" << root4->mActiveVoxelCount << "\n";

    //-------------------------------------------------------------------
    // Hash [x, y, z] to a key.
    // Do a linear search through the root's tiles to find the key (binary search is disabled in nanovdb code).
    //-------------------------------------------------------------------
    struct RootTileData4 {
        uint64_t key;                   // hash key to match
        float    value;                 // value of tile (i.e. no child node)
        int32_t  childID;               // negative values indicate no child node, i.e. this is a value tile
        uint8_t  state;                 // when childID < 0, indicates whether tile is active (1 or 0)
        uint8_t  padding[12];           // to match NanoVDB
    };

    struct RootTileData12 {
        uint64_t key;                   // hash key to match
        float    value[3];              // value of tile (i.e. no child node)
        int32_t  childID;               // negative values indicate no child node, i.e. this is a value tile
        uint8_t  state;                 // when childID < 0, indicates whether tile is active (1 or 0)
        uint8_t  padding[4];            // to match NanoVDB
    };

    const uint root_tile_size = (value_size == 4) ? sizeof( RootTileData4 ) : sizeof( RootTileData12 );

    union InternalTileData4 {
        float    value;                 // value of tile (i.e. no child node)
        int32_t  childID;               // negative values indicate no child node, i.e. this is a value tile
    };

    union InternalTileData12 {
        float    value[3];              // value of tile (i.e. no child node)
        int32_t  childID;               // negative values indicate no child node, i.e. this is a value tile
    };

    const uint8_t * tile_data = root_tile_data;

    die_assert( (32 - INTERNAL2_TOTAL) < 21, "key won't fit in 64 bits" );
    uint64_t key = (uint64_t(uint32_t(z) >> INTERNAL2_TOTAL) <<  0) |
                   (uint64_t(uint32_t(y) >> INTERNAL2_TOTAL) << 21) |
                   (uint64_t(uint32_t(x) >> INTERNAL2_TOTAL) << 42);
    uint i;
    mdout << "VOLUME: root_tile_size=" << root_tile_size << "\n";
    for( i = 0; i < root4->mTileCount; i++ )
    {
        const uint64_t * tile_key = reinterpret_cast< const uint64_t * >( tile_data );
        if ( key == *tile_key ) {
            mdout << "VOLUME: root key=0x" << std::hex << key << std::dec << " tile_i=" << i << "\n";
            break;
        }

        tile_data += root_tile_size;
    }

    if ( i == root4->mTileCount ) {
        // no matching tile => bail out
        mdout << "VOLUME: root key=0x" << std::hex << key << std::dec << " NO_MATCH\n";
        if ( is_active != nullptr ) *is_active = false;
        return nullptr;
    }

    //-------------------------------------------------------------------
    // Check to see if this is a value tile.
    //-------------------------------------------------------------------
    const RootTileData4  * tile4  = reinterpret_cast< const RootTileData4 * >( tile_data );
    const RootTileData12 * tile12 = reinterpret_cast< const RootTileData12 * >( tile_data );
    int32_t childID = (value_size == 4) ? (tile4 ? tile4->childID : -2) : (tile12 ? tile12->childID : -2);
    const uint8_t * value_ptr;
    mdout << "VOLUME: root=" << root4 << " tile=" << tile4 << " childID=" << childID << "\n";
    if ( childID < 0 ) {
        //-------------------------------------------------------------------
        // The value is the same for all voxels in the tile.
        // Point to that.
        //-------------------------------------------------------------------
        value_ptr = tile_data + sizeof( uint64_t );
        if ( is_active != nullptr ) *is_active = tile4->state;

    } else {
        //-------------------------------------------------------------------
        // Get to InternalData for root's childID, which comes after root's tiles.
        //-------------------------------------------------------------------
        const uint internal_tile_size = (value_size == 4) ? sizeof( InternalTileData4 ) : sizeof( InternalTileData12 );
        mdout << "VOLUME: internal_tile_size=" << internal_tile_size << "\n";

        struct InternalData4
        {
            // skip through variable-length stuff
//          MaskT    mValueMask;            
//          MaskT    mChildMask;                   
//          InternalTileData mTable[INTERNAL1_SIZE];
            float    mValueMin; 
            float    mValueMax;
            int32_t  mBBoxMin[3];           // min corner or index bbox
            int32_t  mBBoxMax[3];           // max corner or index bbox
            int32_t  mOffset;               // number of node offsets till first child node
            uint32_t mFlags;                // get this for free (unused)
            uint32_t mPadding[6];           // not sure why this is needed
        };

        struct InternalData12
        {
            // skip through variable-length stuff
//          MaskT    mValueMask;            
//          MaskT    mChildMask;                   
//          InternalTileData mTable[INTERNAL1_SIZE];
            float    mValueMin[3]; 
            float    mValueMax[3];
            int32_t  mBBoxMin[3];           // min corner or index bbox
            int32_t  mBBoxMax[3];           // max corner or index bbox
            int32_t  mOffset;               // number of node offsets till first child node
            uint32_t mFlags;                // get this for free (unused)
            uint32_t mPadding[2];           // not sure why this is needed
        };

        const uint      internal_data_size = (value_size == 4) ? sizeof(InternalData4) : sizeof(InternalData12);
        const uint      internal2_size = 2*INTERNAL2_MASK_SIZE + INTERNAL2_SIZE*internal_tile_size + internal_data_size;
        const uint8_t * internal_data = root_tile_data + root4->mTileCount*root_tile_size + childID*internal2_size;
        const uint8_t * mValueMask = internal_data;
        const uint8_t * mChildMask = mValueMask + INTERNAL2_MASK_SIZE;
        const uint8_t * mTable     = mChildMask + INTERNAL2_MASK_SIZE;
        const InternalData4  * internal4  = reinterpret_cast< const InternalData4 * >(  internal_data + internal2_size - sizeof(InternalData4)  );
        const InternalData12 * internal12 = reinterpret_cast< const InternalData12 * >( internal_data + internal2_size - sizeof(InternalData12) );
              int32_t   mOffset = (value_size == 4) ? internal4->mOffset : internal12->mOffset;
        mdout << "VOLUME: internal_data_size=" << internal_data_size << " internal2_size=" << internal2_size << " LOG2DIM=" << INTERNAL2_LOG2DIM << 
                " mask_size=" << INTERNAL2_MASK_SIZE << " tile_size=" << internal_tile_size << " mtable_size=" << (INTERNAL2_SIZE*internal_tile_size) <<
                " size=" << internal2_size << " mOffset=" << mOffset << "\n";

        //-------------------------------------------------------------------
        // Hash [x,y,z] to child c for internal node.
        // Test mChildMask[c].
        //-------------------------------------------------------------------
        uint64_t c = (((uint64_t(x) & INTERNAL2_MASK) >> INTERNAL1_TOTAL) << (2*INTERNAL2_LOG2DIM)) + 
                     (((uint64_t(y) & INTERNAL2_MASK) >> INTERNAL1_TOTAL) <<    INTERNAL2_LOG2DIM)  +
                     (((uint64_t(z) & INTERNAL2_MASK) >> INTERNAL1_TOTAL) <<                    0);
        tile_data  = mTable + c*internal_tile_size;
        auto tile4   = reinterpret_cast< const InternalTileData4 * >( tile_data );
        auto tile12  = reinterpret_cast< const InternalTileData12 * >( tile_data );
        childID = (value_size == 4) ? (tile4 ? tile4->childID : -2) : (tile12 ? tile12->childID : -2);
        bool is_on = (mChildMask[c/8] >> (c & 7)) & 1;
        mdout << "VOLUME: internal2=" << reinterpret_cast<const void *>(internal_data) << " &mask=" << reinterpret_cast<const void *>(mChildMask) << 
                " n=" << c << " tile=" << tile4 << " childID=" << childID << 
                " isOn=" << is_on << " TOTAL=" << INTERNAL2_TOTAL << " SIZE=" << INTERNAL2_SIZE << "\n";
        if ( !is_on || childID < 0 ) {
            //-------------------------------------------------------------------
            // Pull value out of the mTable[c] Tile.
            //-------------------------------------------------------------------
            value_ptr = tile_data;
            mdout << "VOLUME: &mTable[n]=" << reinterpret_cast<const void *>(tile_data) << 
                           " &mTable[n].value=" << reinterpret_cast<const void *>(value_ptr) << "\n";
            if ( is_active != nullptr ) {
                *is_active = (mValueMask[c/8] >> (c & 7)) & 1; 
                mdout << "VOLUME: is_active=" << *is_active << "\n";
            }

        } else {
            //-------------------------------------------------------------------
            // Recurse to next level of internal nodes.
            // There are only two levels of internal nodes, so we hardcode that for now.
            // If it changes (doubtful), we can make this a loop.
            //
            // Hash [x,y,z] to child c for internal node.
            // Test mChildMask[c].
            //-------------------------------------------------------------------
            const uint internal1_size = 2*INTERNAL1_MASK_SIZE + INTERNAL1_SIZE*internal_tile_size + internal_data_size;
            internal_data += mOffset*internal2_size + childID*internal1_size;
            mValueMask = internal_data;
            mChildMask = mValueMask + INTERNAL1_MASK_SIZE;
            mTable     = mChildMask + INTERNAL1_MASK_SIZE;
            internal4  = reinterpret_cast< const InternalData4 * >(  internal_data + internal1_size - internal_data_size );
            internal12 = reinterpret_cast< const InternalData12 * >( internal_data + internal1_size - internal_data_size );
            mOffset    = (value_size == 4) ? internal4->mOffset : internal12->mOffset;
            mdout << "VOLUME: internal1_size=" << internal1_size << " LOG2DIM=" << INTERNAL1_LOG2DIM << 
                    " mask_size=" << INTERNAL1_MASK_SIZE << " tile_size=" << internal_tile_size << " mtable_size=" << (INTERNAL1_SIZE*internal_tile_size) <<
                    " size=" << internal1_size << " mOffset=" << mOffset << "\n";

            c = (((x & INTERNAL1_MASK) >> LEAF_TOTAL) << (2*INTERNAL1_LOG2DIM)) + 
                (((y & INTERNAL1_MASK) >> LEAF_TOTAL) <<    INTERNAL1_LOG2DIM)  +
                (((z & INTERNAL1_MASK) >> LEAF_TOTAL) <<                    0);
            tile_data  = mTable + c*internal_tile_size;
            tile4  = reinterpret_cast< const InternalTileData4 * >( tile_data );
            tile12 = reinterpret_cast< const InternalTileData12 * >( tile_data );
            childID = (value_size == 4) ? (tile4 ? tile4->childID : -2) : (tile12 ? tile12->childID : -2);
            bool is_on = (mChildMask[c/8] >> (c & 7)) & 1;
            mdout << "VOLUME: internal1=" << reinterpret_cast<const void *>(internal_data) << " &mask=" << reinterpret_cast<const void *>( mChildMask ) << 
                    " n=" << c << " tile=" << reinterpret_cast<const void *>(tile_data) << " childID=" << childID << 
                    " isOn=" << is_on << " TOTAL=" << INTERNAL1_TOTAL << " SIZE=" << INTERNAL1_SIZE << "\n";
            if ( !is_on || childID < 0 ) {
                //-------------------------------------------------------------------
                // Pull value out of the mTable[c] Tile.
                //-------------------------------------------------------------------
                value_ptr = tile_data;
                mdout << "VOLUME: &mTable[n]=" << reinterpret_cast<const void *>(tile_data) << 
                               " &mTable[n].value=" << reinterpret_cast<const void *>(value_ptr) << "\n";
                if ( is_active != nullptr ) {
                    *is_active = (mValueMask[c/8] >> (c & 7)) & 1; 
                    mdout << "VOLUME: is_active=" << *is_active << "\n";
                }

            } else {
                //-------------------------------------------------------------------
                // Recurse to leaf node.
                // Hash [x,y,z] to value index c for leaf node.
                // Pull out mValues[c].
                //-------------------------------------------------------------------
                struct LeafData4
                {
                    // skip through variable-length stuff
//                  MaskT<LEAF_LOG2DIM> mValueMask;      
                    float          mValues[LEAF_SIZE];
                    float          mValueMin;
                    float          mValueMax;
                    _int           mBBoxMin[3]; 
                    uint8_t        mBBoxDif[3]; 
                    uint8_t        mFlags;          // get this for free (unused)
                    uint32_t       mPadding[2];
                };

                struct LeafData12
                {
                    // skip through variable-length stuff
//                  MaskT<LEAF_LOG2DIM> mValueMask;      
                    float          mValues[LEAF_SIZE][3];
                    float          mValueMin[3];
                    float          mValueMax[3];
                    _int           mBBoxMin[3]; 
                    uint8_t        mBBoxDif[3]; 
                    uint8_t        mFlags;          // get this for free (unused)
                    uint32_t       mPadding[6];
                };

                const uint         leaf_data_size = (value_size == 4) ? sizeof(LeafData4) : sizeof(LeafData12);
                const uint         leaf_size  = LEAF_MASK_SIZE + leaf_data_size;
                const uint8_t *    leaf_data  = internal_data + mOffset*internal1_size + childID*leaf_size;
                const uint8_t *    mValueMask = leaf_data;
                const LeafData4 *  leaf4      = reinterpret_cast< const LeafData4 * >(  leaf_data + leaf_size - sizeof(LeafData4) );
                const LeafData12 * leaf12     = reinterpret_cast< const LeafData12 * >( leaf_data + leaf_size - sizeof(LeafData12) );

                c = ((uint64_t(x) & LEAF_MASK) << (2*LEAF_LOG2DIM)) + 
                    ((uint64_t(y) & LEAF_MASK) <<    LEAF_LOG2DIM)  +
                     (uint64_t(z) & LEAF_MASK);
                mdout << "VOLUME: leaf=" << reinterpret_cast<const void *>(leaf_data) << 
                        " &mask=" << reinterpret_cast< const void * >( mValueMask ) << " leaf_size=" << leaf_size <<
                        " n=" << c << "\n";
                const uint8_t * value4_ptr  = reinterpret_cast< const uint8_t * >( &leaf4->mValues[c] );
                const uint8_t * value12_ptr = reinterpret_cast< const uint8_t * >(  leaf12->mValues[c] );
                value_ptr = (value_size == 4) ? value4_ptr : value12_ptr;
                if ( is_active != nullptr ) {
                    *is_active = (mValueMask[c/8] >> (c & 7)) & 1; 
                    mdout << "VOLUME: is_active=" << *is_active << "\n";
                }
            }
        }
    }

    return value_ptr;
}

uint Model::Volume::voxel_class_grid_i( const Model * model, Model::VolumeVoxelClass c ) const
{
    for( uint i = grid_i; i < (grid_i + grid_cnt); i++ )
    {
        if ( model->volume_grids[i].voxel_class == c ) return i;
    }
    return uint(-1);
}

bool Model::Volume::is_emissive( const Model * model ) const
{
    //--------------------------------------------------------------------------
    // Bit masks have been precomputed if this truly is an emissive volume.
    //--------------------------------------------------------------------------
    const uint64_t * super_voxels = model->volume_to_emissive[this - model->volumes];
    return super_voxels != nullptr;
}

bool Model::Volume::rand_emissive_xyz( const Model * model, _int& x, _int& y, _int& z ) const
{
    //--------------------------------------------------------------------------
    // Bit masks have been precomputed if this truly is an emissive volume.
    //--------------------------------------------------------------------------
    const uint64_t * super_voxels = model->volume_to_emissive[this - model->volumes];
    if ( super_voxels == nullptr ) return false;

    //--------------------------------------------------------------------------
    // Start at some random supervoxel and try to find any supervoxel with a non-zero bitmask.
    // If there is such a supervoxel (there should be unless all temperatures are 0 :-), then 
    // pick a random voxel in the supervoxel with its bit set and use that voxel's xyz.
    //--------------------------------------------------------------------------
    uint d_grid_i = voxel_class_grid_i( model, Model::VolumeVoxelClass::DENSITY );
    const VolumeGrid * d_grid_ptr = &model->volume_grids[d_grid_i];
    _int x_min = d_grid_ptr->index_box._min[0];
    _int y_min = d_grid_ptr->index_box._min[0];
    _int z_min = d_grid_ptr->index_box._min[0];
    _int x_max = d_grid_ptr->index_box._max[0];
    _int y_max = d_grid_ptr->index_box._max[0];
    _int z_max = d_grid_ptr->index_box._max[0];
    uint x_cnt = x_max - x_min + 1;
    uint y_cnt = y_max - y_min + 1;
    uint z_cnt = z_max - z_min + 1;
    uint super_x_cnt = (x_cnt + 3) / 4;
    uint super_y_cnt = (y_cnt + 3) / 4;
    uint super_z_cnt = (z_cnt + 3) / 4;
    uint super_cnt = super_x_cnt * super_y_cnt * super_z_cnt;
    uint super_i_start = uint( MODEL_UNIFORM_FN() * float(super_cnt) );
    uint super_i = super_i_start;
    for( ;; )
    {
        if ( super_voxels[super_i] != 0 ) {
            // excellent, we found one
            // pick voxel out of 64 possible
            uint64_t super_voxel = super_voxels[super_i];
            uint b = uint( MODEL_UNIFORM_FN() * 64.0 );
            for( ;; )
            {
                if ( (super_voxel >> b) & 1 ) {
                    // got it
                    z = b / 16;
                    b -= z*16;
                    y = b / 4;
                    x = b - y*4;
                    return true;
                }
                
                // we should eventually find a bit set
                b++;
                if ( b == 16 ) b = 0;  
            }
        }

        // try next supervoxel
        super_i++;
        if ( super_i == super_cnt ) super_i = 0;
        if ( super_i == super_i_start ) return false; // rare
    }
}

inline bool Model::VolumeGrid::is_active( const Model * model, _int x, _int y, _int z ) const
{
    //--------------------------------------------------------------------------
    // Use value_ptr() which can fill in the boolean.
    //--------------------------------------------------------------------------
    bool _is_active;
    value_ptr( model, x, y, z, &_is_active );
    return _is_active;
}

Model::_int64 Model::VolumeGrid::int64_value( const Model * model, _int x, _int y, _int z ) const
{
    auto ptr  = value_ptr( model, x, y, z );
    if ( ptr == nullptr ) return 0;

    _int64 r;
    switch( voxel_type ) 
    {
        case VolumeVoxelType::INT16:    r = *reinterpret_cast<const int16_t *>( ptr ); break;
        case VolumeVoxelType::INT32:    r = *reinterpret_cast<const int32_t *>( ptr ); break;
        case VolumeVoxelType::INT64:    r = *reinterpret_cast<const int64_t *>( ptr ); break;
        default:                        r = 0; die_assert( false, "voxel_type is not INT{16,32,64}" ); break;
    }        
    return r;
}

Model::_int Model::VolumeGrid::int_value( const Model * model, _int x, _int y, _int z ) const
{
    return int64_value( model, x, y, z );
}

Model::real64 Model::VolumeGrid::real64_value( const Model * model, _int x, _int y, _int z ) const
{
    auto ptr  = value_ptr( model, x, y, z );
    if ( ptr == nullptr ) return 0.0;

    real64 r;
    switch( voxel_type ) 
    {
        case VolumeVoxelType::FLOAT:    r = *reinterpret_cast<const float *>( ptr ); break;
        case VolumeVoxelType::DOUBLE:   r = *reinterpret_cast<const double *>( ptr ); break;
        default:                        r = 0; die_assert( false, "voxel_type is not FLOAT or DOUBLE" ); break;
    }        
    return r;
}

Model::real Model::VolumeGrid::real_value( const Model * model, _int x, _int y, _int z ) const
{
    return real64_value( model, x, y, z );
}

Model::real3d Model::VolumeGrid::real3d_value( const Model * model, _int x, _int y, _int z ) const
{
    auto ptr  = value_ptr( model, x, y, z );
    if ( ptr == nullptr ) return real3d( 0, 0, 0 );

    real3d r3d;
    for( int i = 0; i < 3; i++ )
    {
        switch( voxel_type ) 
        {
            case VolumeVoxelType::VEC3F:    r3d.c[i] = *(reinterpret_cast<const float *>( ptr ) + i); break;
            case VolumeVoxelType::VEC3D:    r3d.c[i] = *(reinterpret_cast<const double *>( ptr ) + i); break;
            default:                        r3d.c[i] = 0; die_assert( false, "voxel_type is not VEC3F or VEC3D" ); break;
        }        
    }
    return r3d;
}

Model::real3 Model::VolumeGrid::real3_value( const Model * model, _int x, _int y, _int z ) const
{
    real3d r3d = real3d_value( model, x, y, z );
    real3  r3;
    r3.c[0] = r3d.c[0];
    r3.c[1] = r3d.c[1];
    r3.c[2] = r3d.c[2];
    return r3;
}

bool Model::Volume::hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv,
                         real /*solid_angle*/, real tmin, real tmax, HitRecord& rec ) const
{
    //-------------------------------------------------------------------
    // Clip ray (origin -> direction with t_min,t_max) to our bounding box.
    // Calling AABB::hit() will update tmin,tmax if there's intersection.
    // Then we can get the start and endpoints trivially.
    //-------------------------------------------------------------------
    const VolumeGrid * grid = &model->volume_grids[grid_i];  // all grids in volume have the same bbox
    if ( !grid->world_box.hit( origin, direction, direction_inv, tmin, tmax ) ) { 
        mdout << "Model::Volume::hit: ray does not intersect world_box\n";
        return false;
    }
    if ( tmin < 0.0 )  tmin = 0.0;
    if ( tmax < tmin ) tmax = tmin;

    //-------------------------------------------------------------------
    // Figure out which grids we have.
    // This will affect our decision to stop.
    //-------------------------------------------------------------------
    uint grid_density_i          = -1;
    uint grid_normal_i           = -1;
    uint grid_IOR_i              = -1;
    uint grid_F0_i               = -1;         
    for( uint i = grid_i; i < (grid_i+grid_cnt); i++ )
    {
        switch( model->volume_grids[i].voxel_class ) 
        {
            case VolumeVoxelClass::DENSITY:             grid_density_i = i;             break;
            case VolumeVoxelClass::NORMAL:              grid_normal_i = i;              break;
            case VolumeVoxelClass::IOR:                 grid_IOR_i = i;                 break;
            case VolumeVoxelClass::F0:                  grid_F0_i = i;                  break;
            default:                                                                    break;
        }
    }

    //-------------------------------------------------------------------
    // Start at the clipped ray's origin and voxel.
    // Find a voxel along the ray that has a non-zero value,
    // which we currently interpret like alpha.
    // We use this "alpha" as the probability of stopping at that voxel.
    //-------------------------------------------------------------------
    real64 max_length_within_voxel_inv = 1.0 / grid->world_voxel_size.length();
    mdout << "Model::Volume::hit: model=" << model << " origin=" << origin << " direction=" << direction << 
                                " grid_i=" << grid_i << " grid_density_i=" << grid_density_i << " grid_normal_i=" << grid_normal_i << 
                                " grid_IOR_i=" << grid_IOR_i << " grid_F0_i=" << grid_F0_i << 
                                " world_box=" << grid->world_box << " world_voxel_size=" << grid->world_voxel_size << " world_voxel_size_inv=" << grid->world_voxel_size_inv <<
                                " max_length_within_voxel_inv=" << max_length_within_voxel_inv <<
                                " index_box=" << grid->index_box << " tmin=" << tmin << " tmax=" << tmax << "\n";
    rec.grid_i = grid_i;   // for now
    _int x;
    _int y;
    _int z;
    real voxel_tmin;
    real voxel_tmax;
    real3 p = origin + tmin*direction;
    _int x_prev = -0x7fffffff;
    _int y_prev = -0x7fffffff;
    _int z_prev = -0x7fffffff;
    real t_epsilon = 1e-8;
    for( ;; ) 
    {
        //-------------------------------------------------------------------
        // Get voxel indexes xyz for p.
        // Figure out starting and ending points of ray within the voxel.
        // Calculate length of that line segment.
        //-------------------------------------------------------------------
        real3d pd = real3d( p.c[0], p.c[1], p.c[2] );
        real3d p_scaledd = (pd - grid->world_translate) * grid->world_voxel_size_inv;
        real3  p_scaled = real3( p_scaledd.c[0], p_scaledd.c[1], p_scaledd.c[2] );
        x = std::floor( p_scaled.c[0] );
        y = std::floor( p_scaled.c[1] );
        z = std::floor( p_scaled.c[2] );
        if ( !grid->index_box.encloses( x, y, z ) ) {
            mdout << "Model::Volume::hit: stepped out of volume's index_box=" << grid->index_box << " p_scaled=" << p_scaled << " returning false\n";
            return false;
        }
        real3d voxel_mind = real3d( x, y, z ) * grid->world_voxel_size + grid->world_translate;
        real3  voxel_min  = real3( voxel_mind.c[0], voxel_mind.c[1], voxel_mind.c[2] );
        real3d voxel_maxd = voxel_mind + grid->world_voxel_size;
        real3  voxel_max  = real3( voxel_maxd.c[0], voxel_maxd.c[1], voxel_maxd.c[2] );
        AABB   voxel_box  = AABB( voxel_min );
        voxel_box.expand( voxel_max );
        mdout << "Model::Volume::hit: p=" << p << " p_scaled=" << p_scaled << " xyz=[" << x << "," << y << "," << z << "] voxel_box=" << voxel_box << "\n";
        voxel_tmin = tmin;
        voxel_tmax = tmax;
        if ( (x != x_prev || y != y_prev || z != z_prev) && voxel_box.hit( origin, direction, direction_inv, voxel_tmin, voxel_tmax ) ) { // allowing a little slop here
            if ( voxel_tmin < tmin ) voxel_tmin = tmin;
            if ( voxel_tmax > tmax ) voxel_tmax = tmax;
            real3 voxel_start = origin + voxel_tmin*direction;
            real3 voxel_end   = origin + voxel_tmax*direction;
            real ray_length_within_voxel = (voxel_end-voxel_start).length();
            mdout << "Model::Volume::hit: voxel_tmin=" << voxel_tmin << " voxel_tmax=" << voxel_tmax << 
                                        " voxel_start=" << voxel_start << " voxel_end=" << voxel_end <<
                                        " ray_length_within_voxel=" << ray_length_within_voxel << "\n";

            //-------------------------------------------------------------------
            // See if we should stop.  This occurs if any of these is true:
            //
            // 1) There is no density grid (rare).
            // 2) We have a change in IOR/F0 or if the caller didn't pass in a current F0.
            // 3) We probabilistically choose using the density as the probability, modulated
            //    by mean free path length (MFPL) and length of ray within the voxel.
            //
            // t is just (p - origin) * direction_inv.  We use the component that has the largest value in direction.
            //-------------------------------------------------------------------
            if ( grid_density_i == uint(-1) ) {
                 // treat it like density = 1.0
                 mdout << "Model::Volume::hit: p=" << p << " xyz=[" << x << "," << y << "," << z << "] NO DENSITY\n"; 
                 break;
            }

            if ( grid_F0_i != uint(-1) || grid_IOR_i != uint(-1) ) {
                real3 F0;
                if ( grid_F0_i != uint(-1) ) {
                    rec.grid_i = grid_F0_i;
                    grid = &model->volume_grids[grid_F0_i];
                    F0 = grid->real3_value( model, x, y, z );
                } else {
                    rec.grid_i = grid_IOR_i;
                    grid = &model->volume_grids[grid_IOR_i];
                    real IOR = grid->real_value( model, x, y, z );
                    real a = ((IOR-1.0)*(IOR-1.0)) / ((IOR+1.0)*(IOR+1.0));
                    F0 = real3( a, a, a );
                }
                if ( rec.F0.c[0] < 0.0 || rec.F0.c[0] != F0.c[0] || rec.F0.c[1] != F0.c[1] || rec.F0.c[2] != F0.c[2] ) {
                    mdout << "Model::Volume::hit: p=" << p << " xyz=[" << x << "," << y << "," << z << "] F0 changed to " << F0 << "\n";
                    break;
                }
            }

            grid = &model->volume_grids[grid_density_i];
            real density = grid->real_value( model, x, y, z );
            mdout << "Model::Volume::hit: density=" << density << "\n";
            if ( density > 0.0 ) {
                uint mfpl_grid_i = voxel_class_grid_i(model, Model::VolumeVoxelClass::MFPL);
                real mfpl = (mfpl_grid_i != uint(-1)) ? model->volume_grids[mfpl_grid_i].real_value(model, x, y, z) : 0.1;
                mdout << "Model::Volume::hit: mfpl=" << mfpl << "\n";
                if ( mfpl > 0.0 ) {
                    if ( MODEL_UNIFORM_FN() < (ray_length_within_voxel*max_length_within_voxel_inv*mfpl*density) ) {
                        rec.grid_i = grid_density_i;
                        break;
                    }
                }
            }

            //-------------------------------------------------------------------
            // Keep going just within the next voxel.
            //-------------------------------------------------------------------
            t_epsilon = (voxel_tmax-voxel_tmin)*0.1;  
            if ( t_epsilon <= 0.0 ) t_epsilon = 1e-8;
        } else {
            //-------------------------------------------------------------------
            // Double t_epsilon and try again.
            // We are probably very close to a corner.
            //-------------------------------------------------------------------
            t_epsilon *= 2.0;
        }

        real t_new = voxel_tmax + t_epsilon;
        mdout << "Model::Volume::hit: t_epsilon=" << t_epsilon << " t_new=" << t_new << "\n";
        if ( t_new > tmax ) {
            mdout << "Model::Volume::hit: t_epsilon move stepped out of volume, returning false\n";
            return false;
        }
        p = origin + t_new*direction;
        x_prev = x;
        y_prev = y;
        z_prev = z;
    } // for

    //-------------------------------------------------------------------
    // Change p to some uniformly random location between voxel_start and voxel_end.
    //-------------------------------------------------------------------
    real t = voxel_tmin + MODEL_UNIFORM_FN()*(voxel_tmax-voxel_tmin);
    p = origin + t*direction;
    
    //-------------------------------------------------------------------
    // Fill in rest of rec.
    //-------------------------------------------------------------------
    rec.model = model;
    rec.p = real3( p.c[0], p.c[1], p.c[2] );
    int c;
    for( int i=0; i < 3; i++ )
    {
        if ( i == 0 || direction.c[i] > direction.c[c] ) {
            c = i;
        }
    }
    rec.t = (p.c[c] - origin.c[c]) * direction_inv.c[c];
    rec.poly_i = uint(-1);
    rec.volume_i = this - model->volumes;
    rec.voxel_xyz[0] = x;
    rec.voxel_xyz[1] = y;
    rec.voxel_xyz[2] = z;
    grid = &model->volume_grids[rec.grid_i];
    switch( grid->voxel_type )
    {
        case VolumeVoxelType::INT16:            rec.voxel_value.i   = grid->int_value( model, x, y, z );           break;
        case VolumeVoxelType::INT32:            rec.voxel_value.i   = grid->int_value( model, x, y, z );           break;
        case VolumeVoxelType::INT64:            rec.voxel_value.i64 = grid->int64_value( model, x, y, z );         break;
        case VolumeVoxelType::FP16:             rec.voxel_value.r   = grid->real_value( model, x, y, z );          break;
        case VolumeVoxelType::FLOAT:            rec.voxel_value.r   = grid->real_value( model, x, y, z );          break;
        case VolumeVoxelType::DOUBLE:           rec.voxel_value.r64 = grid->real64_value( model, x, y, z );        break;
        case VolumeVoxelType::VEC3F:            rec.voxel_value.r3  = grid->real3_value( model, x, y, z );         break;
        case VolumeVoxelType::VEC3D:            rec.voxel_value.r3d = grid->real3d_value( model, x, y, z );        break;
        default:                                                                                                        break;
    }
    if ( grid_normal_i != uint(-1) ) {
        rec.normal = model->volume_grids[grid_normal_i].real3_value( model, x, y, z );
    } else {
        rec.normal = real3(0,1,0);
    }
    rec.shading_normal = rec.normal;
    mdout << "Model::Volume::hit: success origin=" << origin << " direction=" << direction << " direction_inv=" << direction_inv << 
             " p=" << p << " c=" << c << " t=" << rec.t << "\n";
    return true;  // success
}

bool Model::Instance::bounding_box( const Model * model, AABB& b, real padding ) const
{
    (void)model;
    b = box;
    b.pad( padding );
    return true;
}

bool Model::Instance::hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                           real solid_angle, real t_min, real t_max, HitRecord& rec )
{
    die_assert( sizeof(real) == 4, "real is not float" );
    die_assert( kind == INSTANCE_KIND::MODEL_PTR, "did not get model_ptr instance" );

    //-----------------------------------------------------------
    // If the instance is empty, it can have a null BVH.
    //-----------------------------------------------------------
    Model* t_model = u.model_ptr;
    if ( t_model->hdr->bvh_root_i == uint(-1) ) return false;

    //-----------------------------------------------------------
    // Use inverse matrix to transform origin, direction, and direction_inv into target model's space.
    // Then call model's root BVH hit node.
    //-----------------------------------------------------------
    Matrix * M_inv = &model->matrixes[matrix_inv_i];
    real3    t_origin, t_direction, t_direction_inv;
    M_inv->transform( origin, t_origin );
    M_inv->transform( direction, t_direction );
    M_inv->transform( direction_inv, t_direction_inv );
    mdout << "Model::Instance::hit: model=" << model << " M_inv=" << *M_inv << 
                                  " origin="   << origin   << " direction="   << direction   << " direction_inv="   << direction_inv << 
                                  " t_model=" << t_model << " t_origin=" << t_origin << " t_direction=" << t_direction << " t_direction_inv=" << t_direction_inv << "\n";

    BVH_Node * t_bvh = &t_model->bvh_nodes[t_model->hdr->bvh_root_i];
    if ( !t_bvh->hit( t_model, t_origin, t_direction, t_direction_inv, solid_angle, t_min, t_max, rec ) ) return false;

    //-----------------------------------------------------------
    // Use matrix to transform rec.p back to global world space.
    // Use transposed inverse matrix to transform rec.normal correctly (ask Pete Shirley).
    //-----------------------------------------------------------
    Matrix * M           = &model->matrixes[matrix_i];
    Matrix * M_inv_trans = &model->matrixes[matrix_inv_trans_i];
    real3 p = rec.p;
    real3 normal = rec.normal;
    M->transform( p, rec.p );
    M_inv_trans->transform( normal, rec.normal );
    return true;
}

inline bool Model::BVH_Node::bounding_box( const Model * model, Model::AABB& b ) const
{
    (void)model;
    b = box;
    return true;
}

void Model::bvh_build( Model::BVH_TREE bvh_tree )
{
    //--------------------------------------------------
    // Clear out any prior BVH.
    //--------------------------------------------------
    hdr->bvh_node_cnt = 0;
    hdr->bvh_root_i = uint(-1);
    if ( bvh_tree == BVH_TREE::NONE ) return;

    if ( hdr->poly_cnt != 0 ) {
        //--------------------------------------------------
        // Build a proper BVH with polygons at the leaves
        //--------------------------------------------------
        die_assert( hdr->inst_cnt == 0, "inst_cnt should be 0 when poly_cnt is not 0" );
        hdr->bvh_root_i = bvh_node( BVH_NODE_KIND::POLYGON, 0, hdr->poly_cnt, 1 );

        if ( hdr->emissive_poly_cnt != 0 ) {
            // need to rebuild this due to changes in poly indexes
            // should end up with the same number of emissive polys
            //
            uint epoly_cnt = 0;
            for( uint i = 0; i < hdr->poly_cnt; i++ )
            {
                uint mtl_i = polygons[i].mtl_i;
                if ( mtl_i != uint(-1) && materials[mtl_i].is_emissive() ) {
                    emissive_polygons[epoly_cnt++] = i;
                }
            }
            if ( epoly_cnt != hdr->emissive_poly_cnt ) {
                std::cout << "ERROR: after BVH build, new emissive_poly_cnt=" << epoly_cnt << 
                             " != old emissive_poly_cnt=" << hdr->emissive_poly_cnt << "\n";
                exit( 1 );
            }
        }
    } else if ( hdr->volume_cnt != 0 ) {
        //--------------------------------------------------
        // For grids.
        //--------------------------------------------------
        die_assert( hdr->inst_cnt == 0, "inst_cnt should be 0 when volume_cnt is not 0" );
        hdr->bvh_root_i = bvh_node( BVH_NODE_KIND::VOLUME, 0, hdr->volume_cnt, 1 );

    } else if ( hdr->inst_cnt != 0 ) {
        //--------------------------------------------------
        // For instances.
        //--------------------------------------------------
        hdr->bvh_root_i = bvh_node( BVH_NODE_KIND::INSTANCE, 0, hdr->inst_cnt, 1 );
    }
}

inline Model::uint Model::bvh_qsplit( BVH_NODE_KIND kind, Model::uint first, Model::uint n, Model::real pivot, Model::uint axis )
{
    uint m = first;

    for( uint i = first; i < (first+n); i++ )
    {
        AABB box;
        polygons[i].bounding_box( this, box );
        real centroid = (box._min.c[axis] + box._max.c[axis]) * 0.5;
        if ( centroid < pivot ) {
            switch( kind )
            {
                case BVH_NODE_KIND::POLYGON:
                {
                    Polygon temp = polygons[i];
                    polygons[i]  = polygons[m];
                    polygons[m]  = temp;
                    break;
                }
                case BVH_NODE_KIND::VOLUME:
                {
                    Volume temp = volumes[i];
                    volumes[i]  = volumes[m];
                    volumes[m]  = temp;
                    break;
                }
                case BVH_NODE_KIND::INSTANCE:
                {
                    Instance temp = instances[i];
                    instances[i]  = instances[m];
                    instances[m]  = temp;
                    break;
                }
                default:
                {
                    die_assert( false, "unexpected BVH_NODE_KIND" );
                    break;
                }
            }
            m++;
         }
     }

     die_assert( m >= first && m <= (first+n), "qsplit has gone mad" );
     if ( m == first || m == (first +n) ) m = first + n/2;
     return m;
}

Model::uint Model::bvh_node( BVH_NODE_KIND kind, Model::uint first, Model::uint n, Model::uint axis ) 
{
    perhaps_realloc( bvh_nodes, hdr->bvh_node_cnt, max->bvh_node_cnt, 1 );
    uint bvh_i = hdr->bvh_node_cnt++;
    BVH_Node * node = &bvh_nodes[bvh_i];

    switch( kind )
    {
        case BVH_NODE_KIND::POLYGON:                polygons[first+0].bounding_box( this, node->box );          break;
        case BVH_NODE_KIND::VOLUME:                 volumes[first+0].bounding_box( this, node->box );           break;
        case BVH_NODE_KIND::INSTANCE:               instances[first+0].bounding_box( this, node->box );         break;
        default:                                    die_assert( false, "unexpected BVH_NODE_KIND" );            break;
    }

    AABB new_box;
    for( uint i = 1; i < n; i++ )
    {
        switch( kind )
        {
            case BVH_NODE_KIND::POLYGON:                polygons[first+i].bounding_box( this, new_box );        break;
            case BVH_NODE_KIND::VOLUME:                 volumes[first+i].bounding_box( this, new_box );         break;
            case BVH_NODE_KIND::INSTANCE:               instances[first+i].bounding_box( this, new_box );       break;
            default:                                    die_assert( false, "unexpected BVH_NODE_KIND" );        break;
        }
        node->box.expand( new_box );
    }

    if ( n == 1 || n == 2 ) {
        node->left_i = first;
        node->left_kind  = kind;
        node->right_kind = kind;

        switch( kind )
        {
            case BVH_NODE_KIND::POLYGON:                polygons[first+0].bounding_box( this, new_box );        break;
            case BVH_NODE_KIND::VOLUME:                 volumes[first+0].bounding_box( this, new_box );         break;
            case BVH_NODE_KIND::INSTANCE:               instances[first+0].bounding_box( this, new_box );       break;
            default:                                    die_assert( false, "unexpected BVH_NODE_KIND" );        break;
        }
        die_assert( node->box.encloses( new_box ), "box should enclose new_box" );
        if ( n == 2 ) {
            switch( kind )
            {
                case BVH_NODE_KIND::POLYGON:                polygons[first+1].bounding_box( this, new_box );        break;
                case BVH_NODE_KIND::VOLUME:                 volumes[first+1].bounding_box( this, new_box );         break;
                case BVH_NODE_KIND::INSTANCE:               instances[first+1].bounding_box( this, new_box );       break;
                default:                                    die_assert( false, "unexpected BVH_NODE_KIND" );        break;
            }
            die_assert( node->box.encloses( new_box ), "box should enclose new_box" );
            node->right_i = first + 1;
        } else {
            node->right_i = first;
        }

    } else {
        node->left_kind = Model::BVH_NODE_KIND::BVH_NODE;
        node->right_kind = Model::BVH_NODE_KIND::BVH_NODE;
        real pivot = (node->box._min.c[axis] + node->box._max.c[axis]) * 0.5;
        uint m = bvh_qsplit( kind, first, n, pivot, axis );
        uint nm = m - first;
        uint left_i  = bvh_node( kind, first, nm,   (axis + 1) % 3 );
        uint right_i = bvh_node( kind, m,     n-nm, (axis + 1) % 3 );
        node = &bvh_nodes[bvh_i];  // could change after previous calls
        node->left_i  = left_i;
        node->right_i = right_i;
        die_assert( node->box.encloses( bvh_nodes[left_i].box ) && node->box.encloses( bvh_nodes[right_i].box ), "box does not enclose left and right boxes" );
    }

    return bvh_i;
}

inline bool Model::BVH_Node::hit( const Model * model, const Model::real3& origin, 
                                  const Model::real3& direction, const Model::real3& direction_inv,
                                  Model::real solid_angle, Model::real t_min, Model::real t_max, Model::HitRecord& rec ) const
{
    bool r = false;
    uint bvh_i = this - model->bvh_nodes;
    mdout << "Model::BVH_Node::hit: model=" << model << " bvh_i=" << bvh_i << "\n";
    real tmin = t_min;
    real tmax = t_max;
    if ( box.hit( origin, direction, direction_inv, tmin, tmax ) ) {
        HitRecord left_rec;
        HitRecord right_rec;
        bool hit_left  = (left_kind == BVH_NODE_KIND::POLYGON)      ? model->polygons[left_i].hit(      model, origin, direction, direction_inv,   
                                                                                                        solid_angle, t_min, t_max, left_rec  ) :
                         (left_kind == BVH_NODE_KIND::VOLUME)       ? model->volumes[left_i].hit(       model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, left_rec  ) :
                         (left_kind == BVH_NODE_KIND::INSTANCE)     ? model->instances[left_i].hit(     model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, left_rec  ) :
                                                                      model->bvh_nodes[left_i].hit(     model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, left_rec  );
        bool hit_right = (left_i == right_i)                        ? false :   // lone leaf
                         (right_kind == BVH_NODE_KIND::POLYGON)     ? model->polygons[right_i].hit(     model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, right_rec ) :
                         (right_kind == BVH_NODE_KIND::VOLUME)      ? model->volumes[right_i].hit(      model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, right_rec ) :
                         (right_kind == BVH_NODE_KIND::INSTANCE)    ? model->instances[right_i].hit(    model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, right_rec ) :
                                                                      model->bvh_nodes[right_i].hit(    model, origin, direction, direction_inv, 
                                                                                                        solid_angle, t_min, t_max, right_rec );
        if ( hit_left && (!hit_right || left_rec.t < right_rec.t) ) {
            rec = left_rec;
            r = true;
        } else if ( hit_right ) {
            rec = right_rec;
            r = true;
        }
    }
    mdout << "Model::BVH_Node::hit: bvh_i=" << bvh_i << " return=" << r << "\n";
    return r;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// PARSING UTILITIES
//
// These assume an in-memory copy of the entire file with xxx pointing to the 
// current character and xxx_end pointing one character past the last character
// in the file.
//
// For speed, these are self-contained and don't use any built-in C++ functions (slow).
//
// These functions can be used in non-graphics applications.
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline bool Model::skip_whitespace( const char *& xxx, const char * xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx != mtl ) line_num++;
            if ( Model::debug ) {
                std::string line = "";
                for( const char * ccc = xxx+1; ccc != xxx_end && *ccc != '\n' && *ccc != '\r'; ccc++ )
                {
                    line += std::string( 1, *ccc );
                }
                mdout << "PARSING: " << line << "\n";
            }
            in_comment = false;
        }
        xxx++;
    }
    return true;
}

inline bool Model::skip_whitespace_to_eol( const char *& xxx, const char * xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx != mtl ) line_num++;
            break;
        }
        xxx++;
    }
    return true;
}

inline bool Model::skip_to_eol( const char *& xxx, const char * xxx_end )
{
    if ( !eol( xxx, xxx_end ) ) {
        while( xxx != xxx_end )
        {
            char ch = *xxx;
            if ( ch == '\n' || ch == '\r' ) break;
            xxx++;
        }
    }
    return true;
}

inline bool Model::skip_through_eol( const char *& xxx, const char * xxx_end )
{
    while( !eol( xxx, xxx_end ) ) 
    {
        xxx++;
    }
    return true;
}

inline bool Model::eol( const char *& xxx, const char * xxx_end )
{
    skip_whitespace_to_eol( xxx, xxx_end );

    if ( xxx == xxx_end || *xxx == '\n' || *xxx == '\r' ) {
        if ( xxx != xxx_end ) {
            if ( *xxx == '\n' && xxx != mtl ) line_num++;
            xxx++;
        }
        dprint( "at eol" );
        return true;
    } else {
        dprint( "not at eol, char='" + std::string( 1, *xxx ) + "'" );
        return false;
    }
}

inline bool Model::expect_char( char ch, const char *& xxx, const char * xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );
    rtn_assert( xxx != xxx_end, "premature end of file" );
    rtn_assert( *xxx == ch, "expected character '" + std::string(1, ch) + "' got '" + std::string( 1, *xxx ) + "' " + surrounding_lines( xxx, xxx_end ) );
    xxx++;
    return true;
}

inline bool Model::expect_eol( const char *& xxx, const char * xxx_end )
{
    if ( xxx != xxx_end ) {
        rtn_assert( *xxx == '\n' || *xxx == '\r', "not at eol" );
        xxx++;
    }
    return true;
}

inline bool Model::expect_cmd( const char * s, const char *& xxx, const char * xxx_end )
{
    char s_ch1 = '!';
    while( xxx != xxx_end ) 
    {
        s_ch1 = *s;
        s++;
        if ( s_ch1 == '\0' ) {
            // command needs to end with space
            rtn_assert( *xxx == ' ' || *xxx == '\t', "unknown command" );
            return true;
        }

        // case insensitive
        char s_ch2;
        if ( s_ch1 >= 'a' && s_ch1 <= 'z' ) {
            s_ch2 = 'A' + s_ch1 - 'a';
        } else {
            s_ch2 = s_ch1;
        }

        char o_ch = *xxx;
        xxx++;
        rtn_assert( o_ch == s_ch1 || o_ch == s_ch2, "unexpected command character: " + std::string( 1, o_ch ) );
    }

    return s_ch1 == '\0';
}

inline bool Model::parse_string( std::string& s, const char *& xxx, const char * xxx_end )
{
    if ( !expect_char( '"', xxx, xxx_end, true ) ) return false;
    s = "";
    for( ;; ) 
    {
        rtn_assert( xxx != xxx_end, "no terminating \" for string" );
        if ( *xxx == '"' ) {
            xxx++;
            return true;
        }
        s += *xxx;
        xxx++;
    }
}

inline bool Model::parse_string_i( uint& s_i, const char *& xxx, const char * xxx_end )
{
    std::string s;
    if ( !parse_string( s, xxx, xxx_end ) ) return false;

    s_i = make_string( s );
    return true;
}

inline bool Model::parse_name( const char *& name, const char *& xxx, const char * xxx_end )
{
    bool vld = false;
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, 1024 );
    char * _name = &strings[hdr->char_cnt];
    name = _name;

    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    uint len = 0;
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( ch == '\n' || ch == '\r' ) break;

        rtn_assert( len < 1024, "string is larger than 1024 characters" );
        _name[len++] = ch;
        vld = true;
        xxx++;
    }

    if ( vld ) {
        _name[len] = '\0';
        hdr->char_cnt += len+1;
        char * ptr;
        for( ptr = &_name[len-1]; ptr != _name; ptr-- )
        {
            // skip trailing spaces
            if ( *ptr != ' ' && *ptr != '\t' ) break;
            *ptr = '\0';
        }

        return *ptr != '\0';
    }

    rtn_assert( 0, "could not parse name: " + surrounding_lines( xxx, xxx_end ) );
}

inline bool Model::parse_id( std::string& id, const char *& xxx, const char * xxx_end )
{
    skip_whitespace( xxx, xxx_end );

    id = "";
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( !( ch == '_' || (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (id != "" && ch >= '0' && ch <= '9')) ) break;

        id += std::string( 1, ch );
        xxx++;
    }

    return true;
}

inline bool Model::parse_id_i( uint& id_i, const char *& xxx, const char * xxx_end )
{
    std::string id;
    if ( !parse_id( id, xxx, xxx_end ) ) return false;

    uint id_len = id.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, id_len+1 );
    id_i = hdr->char_cnt;
    char * to_id = &strings[hdr->char_cnt];
    hdr->char_cnt += id_len + 1;
    memcpy( to_id, id.c_str(), id_len+1 );
    return true;
}

inline bool Model::parse_option_name( std::string& option_name, const char *& xxx, const char * xxx_end )
{
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces
    if ( xxx == xxx_end or *xxx != '-' ) return false;
    xxx++;

    option_name = "-";
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( !( ch == '_' || (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ) ) break;

        option_name += std::string( 1, ch );
        xxx++;
    }

    return true;
}

inline bool Model::parse_obj_cmd( obj_cmd_t& cmd )
{
    rtn_assert( obj != obj_end, "no .obj command" );

    char ch = *obj;
    obj++;
    switch( ch )
    {
        case 'o':           
        case 'O':           
            rtn_assert( obj != obj_end && *obj == ' ', "bad .obj o command" );
            cmd = CMD_O;
            return true;

        case 'g':           
        case 'G':           
            rtn_assert( obj != obj_end && *obj == ' ', "bad .obj g command" );
            cmd = CMD_G;
            return true;

        case 's':           
        case 'S':           
            rtn_assert( obj != obj_end && *obj == ' ', "bad .obj s command" );
            cmd = CMD_S;
            return true;

        case 'f':           
        case 'F':           
            rtn_assert( obj != obj_end && *obj == ' ', "bad .obj f command" );
            cmd = CMD_F;
            return true;

        case 'v':
        case 'V':
            rtn_assert( obj != obj_end, "truncated .obj v command" );
            ch = *obj;
            obj++;

            if ( ch == 'n' || ch == 'N' ) {
                rtn_assert( obj != obj_end && *obj == ' ', "bad .obj vn command" );
                obj++;
                cmd = CMD_VN;
                return true;
            }
                
            if ( ch == 't' || ch == 'T' ) {
                rtn_assert( obj != obj_end && *obj == ' ', "bad .obj vt command" );
                obj++;
                cmd = CMD_VT;
                return true;
            }

            if ( ch == ' ' ) {
                cmd = CMD_V;
                return true;
            }

            rtn_assert( 0, "bad .obj V command" );

        case 'm':
            cmd = CMD_MTLLIB;
            return expect_cmd( "tllib", obj, obj_end ); 
            
        case 'u':
            cmd = CMD_USEMTL;
            return expect_cmd( "semtl", obj, obj_end ); 
            
        default:
            rtn_assert( 0, "bad .obj command character: " + surrounding_lines( obj, obj_end ) );
    }
}

inline bool Model::parse_mtl_cmd( mtl_cmd_t& cmd, const char *& mtl, const char * mtl_end )
{
    rtn_assert( mtl != mtl_end, "no .mtl command" );

    char ch = *mtl;
    mtl++;
    switch( ch )
    {
        case 'n':
        case 'N':
            rtn_assert( mtl != mtl_end, "truncated .mtl v command" );
            ch = *mtl;
            mtl++;

            if ( ch == 's' || ch == 'S' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Ns command" );
                cmd = CMD_NS;
                return true;
            }

            if ( ch == 'i' || ch == 'I' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Ni command" );
                cmd = CMD_NI;
                return true;
            }

            if ( ch == 'e' || ch == 'E' ) {
                cmd = CMD_NEWMTL;
                return expect_cmd( "wmtl", mtl, mtl_end );
            }

            rtn_assert( 0, "bad .mtl N command" );

        case 'i':
        case 'I':
            cmd = CMD_ILLUM;
            return expect_cmd( "llum", mtl, mtl_end );

        case 'd':           
        case 'D':           
            rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl d command" );
            mtl++;
            cmd = CMD_D;
            return true;

        case 'k':
        case 'K':
            rtn_assert( mtl != mtl_end, "truncated .mtl K command" );
            ch = *mtl;
            mtl++;

            if ( ch == 'a' || ch == 'A' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Ka command" );
                mtl++;
                cmd = CMD_KA;
                return true;
            }
                
            if ( ch == 'd' || ch == 'D' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Kd command" );
                mtl++;
                cmd = CMD_KD;
                return true;
            }

            if ( ch == 'e' || ch == 'E' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Ke command" );
                mtl++;
                cmd = CMD_KE;
                return true;
            }

            if ( ch == 's' || ch == 'S' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Ks command" );
                mtl++;
                cmd = CMD_KS;
                return true;
            }

            rtn_assert( 0, "truncated .mtl K command" );

        case 't':
        case 'T':
            rtn_assert( mtl != mtl_end, "truncated .mtl T command" );
            ch = *mtl;
            mtl++;

            if ( ch == 'f' || ch == 'F' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Tf command" );
                mtl++;
                cmd = CMD_TF;
                return true;
            }
                
            if ( ch == 'r' || ch == 'R' ) {
                rtn_assert( mtl != mtl_end && *mtl == ' ', "bad .mtl Tr command" );
                mtl++;
                cmd = CMD_TR;
                return true;
            }

            rtn_assert( 0, "truncated .mtl T command" );

        case 'm':
        case 'M':
            if ( !expect_char( 'a', mtl, mtl_end ) || 
                 !expect_char( 'p', mtl, mtl_end ) || 
                 !expect_char( '_', mtl, mtl_end ) ) return false;

            rtn_assert( mtl != mtl_end, "truncated .mtl map_ command" );
            ch = *mtl;
            mtl++;

            if ( ch == 'k' || ch == 'K' ) {
                rtn_assert( mtl != mtl_end, "truncated .mtl map_k command" );

                ch = *mtl;
                if ( ch == 'a' || ch == 'A' ) {
                    cmd = CMD_MAP_KA;
                } else if ( ch == 'd' || ch == 'D' ) {
                    cmd = CMD_MAP_KD;
                } else if ( ch == 'e' || ch == 'E' ) { 
                    cmd = CMD_MAP_KE;
                } else if ( ch == 's' || ch == 'S' ) { 
                    cmd = CMD_MAP_KS;
                } else {
                    rtn_assert( 0, "bad .mtl map_k command" );
                }

                mtl++;
                if ( !expect_char( ' ', mtl, mtl_end ) ) {
                    rtn_assert( 0, "bad .mtl map_k command (no space after Ka/Kd/Ke/Ks)" );
                }
                return true;

            } else if ( ch == 'n' || ch == 'N' ) {
                if ( !expect_char( 's', mtl, mtl_end ) || 
                     !expect_char( ' ', mtl, mtl_end ) ) {
                    rtn_assert( 0, "unexpected .mtl map_n command" );
                }
                cmd = CMD_MAP_NS;
                return true;

            } else if ( ch == 'd' || ch == 'D' ) {
                if ( !expect_char( ' ', mtl, mtl_end ) ) {
                    rtn_assert( 0, "unexpected .mtl map_d command" );
                }
                cmd = CMD_MAP_D;
                return true;

            } else if ( ch == 'b' || ch == 'B' ) {
                if ( !expect_char( 'u', mtl, mtl_end ) || 
                     !expect_char( 'm', mtl, mtl_end ) || 
                     !expect_char( 'p', mtl, mtl_end ) ||
                     !expect_char( ' ', mtl, mtl_end ) ) {
                    rtn_assert( 0, "unexpected .mtl map_b command" );
                }

                cmd = CMD_MAP_BUMP;
                return true;

            } else if ( ch == 'r' || ch == 'R' ) {
                if ( !expect_char( 'e', mtl, mtl_end ) || 
                     !expect_char( 'f', mtl, mtl_end ) || 
                     !expect_char( 'l', mtl, mtl_end ) ||
                     !expect_char( ' ', mtl, mtl_end ) ) {
                    rtn_assert( 0, "unexpected .mtl map_b command" );
                }

                cmd = CMD_MAP_REFL;
                return true;
            }

            rtn_assert( 0, "truncated .mtl m command: " + surrounding_lines( mtl, mtl_end ) );

        case 'b':
        case 'B':
            if ( !expect_char( 'u', mtl, mtl_end ) || 
                 !expect_char( 'm', mtl, mtl_end ) || 
                 !expect_char( 'p', mtl, mtl_end ) ||
                 !expect_char( ' ', mtl, mtl_end ) ) {
                rtn_assert( 0, "unexpected .mtl b command" );
            }

            cmd = CMD_MAP_BUMP;
            return true;

        case 'r':
        case 'R':
            if ( !expect_char( 'e', mtl, mtl_end ) || 
                 !expect_char( 'f', mtl, mtl_end ) || 
                 !expect_char( 'l', mtl, mtl_end ) ||
                 !expect_char( ' ', mtl, mtl_end ) ) {
                rtn_assert( 0, "unexpected .mtl b command" );
            }

            cmd = CMD_MAP_REFL;
            return true;

        default:
            rtn_assert( 0, "bad .mtl command character '" + std::string( 1, ch ) + "': " + surrounding_lines( mtl, mtl_end ) );
    }
}

inline bool Model::parse_gph_real( uint& r_i, const char *& gph, const char * gph_end )
{
    real r;
    if ( !parse_real( r, gph, gph_end, true ) ) return false;

    r_i = gph_node_alloc( GRAPH_NODE_KIND::REAL );
    mdout << "    " << r << "\n";
    graph_nodes[r_i].u.r = r;

    return true;
}

inline bool Model::parse_gph_real3( uint& r3_i, const char *& gph, const char * gph_end )
{
    real3 r3;
    if ( !parse_real3( r3, gph, gph_end, true ) ) return false;

    r3_i = gph_node_alloc( GRAPH_NODE_KIND::REAL3 );
    mdout << "    " << r3 << "\n";
    graph_nodes[r3_i].u.r3 = r3;

    return true;
}

inline bool Model::parse_gph_str( uint& s_i, const char *& gph, const char * gph_end )
{
    uint str_i;
    if ( !parse_string_i( str_i, gph, gph_end ) ) return false;

    s_i = gph_node_alloc( GRAPH_NODE_KIND::STR );
    mdout << "    " << &strings[str_i] << "\n";
    graph_nodes[s_i].u.s_i = str_i;

    return true;
}

inline bool Model::parse_gph_id( uint& id_i, const char *& gph, const char * gph_end )
{
    uint s_i;
    if ( !parse_id_i( s_i, gph, gph_end ) ) return false;

    id_i = gph_node_alloc( GRAPH_NODE_KIND::ID );
    mdout << "    " << &strings[s_i] << "\n";
    graph_nodes[id_i].u.s_i = s_i;

    return true;
}

inline bool Model::parse_gph_op_lookahead( uint& s_i, const char *& gph, const char * gph_end )
{
    skip_whitespace( gph, gph_end );
    if ( gph == gph_end ) return false;

    char s[2] = { *gph, '\0' };
    uint s_len = 1;
    switch( s[0] )
    {
        case '=':     
        case '+':    
        case '-':     
        case '*':    
        case '/':   
            break;

        default:
            return false;
    }

    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    s_i = hdr->char_cnt;
    char * to_s = &strings[hdr->char_cnt];
    hdr->char_cnt += s_len+1;
    memcpy( to_s, s, s_len+1 );
    return true;
}

inline bool Model::parse_gph_op( uint& op_i, const char *& gph, const char * gph_end )
{
    uint s_i;
    if ( !parse_gph_op_lookahead( s_i, gph, gph_end ) ) return false;
    gph += strlen( &strings[s_i] );

    op_i = gph_node_alloc( GRAPH_NODE_KIND::OP );
    graph_nodes[op_i].u.s_i = s_i;
    mdout << "    " << std::string(&strings[s_i]) << "\n";

    return true;
}

bool Model::parse_gph_term( uint& term_i, const char *& gph, const char * gph_end )
{
    skip_whitespace( gph, gph_end );
    char ch = *gph;
    if ( ch >= '0' && ch <= '9' ) {
        // constant real
        if ( !parse_gph_real( term_i, gph, gph_end ) ) return false;

    } else if ( ch == '[' ) {
        // constant real3
        if ( !parse_gph_real3( term_i, gph, gph_end ) ) return false;

    } else if ( ch == '"' ) {
        // constant string
        if ( !parse_gph_str( term_i, gph, gph_end ) ) return false;

    } else if ( ch == '(' ) {
        // parenthisized
        gph++;
        if ( !parse_gph_expr( term_i, gph, gph_end ) ) return false;
        skip_whitespace( gph, gph_end );
        rtn_assert( gph != gph_end, "missing ')' at end of .gph file " + surrounding_lines( gph, gph_end ) );
        rtn_assert( *gph == ')', "missing ')' in .gph file , got '" + std::string( gph ) + "' (" + std::to_string(*gph) + ")" +
                                 surrounding_lines( gph, gph_end ) );
        gph++;

    } else if ( (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || ch == '_' ) {
        // identifier
        if ( !parse_gph_id( term_i, gph, gph_end ) ) return false;
        skip_whitespace( gph, gph_end );
        if ( *gph == '(' ) {
            // function call -> parse arguments
            gph++;

            uint id_i = term_i;
            term_i = gph_node_alloc( GRAPH_NODE_KIND::NARY );
            graph_nodes[term_i].u.child_first_i = id_i;
            graph_nodes[id_i].parent_i = term_i;

            uint last_i = id_i;
            uint arg_i;
            for( ;; ) 
            {
                skip_whitespace( gph, gph_end );
                rtn_assert( gph != gph_end, "missing arguments and ')' in function call" + surrounding_lines( gph, gph_end ) );
                if ( *gph == ')' ) {
                    gph++;
                    break;
                }

                if ( last_i != id_i ) {
                    skip_whitespace( gph, gph_end );
                    rtn_assert( *gph == ',', "missing comma in arguments to function call" + surrounding_lines( gph, gph_end ) );
                    gph++;
                }

                if ( !parse_gph_expr( arg_i, gph, gph_end ) ) return false;
                graph_nodes[last_i].sibling_i = arg_i;
                graph_nodes[arg_i].parent_i = term_i;
                last_i = arg_i;
            }
        }
    } else {
        rtn_assert( false, "could not parse gph term, got character '" + std::string( 1, ch ) + surrounding_lines( gph, gph_end ) );
    }

    return true;
}

inline bool Model::parse_gph_expr( uint& expr_i, const char *& gph, const char * gph_end, uint curr_prec )
{
    //-------------------------------------------------------------------
    // Parse one term
    //-------------------------------------------------------------------
    if ( !parse_gph_term( expr_i, gph, gph_end ) ) return false;

    //-------------------------------------------------------------------
    // Now see if we get a bunch of ops in a row that we can parse now.
    //-------------------------------------------------------------------
    for( ;; ) 
    {
        //-------------------------------------------------------------------
        // See if there is an operator next that we can parse now.
        //-------------------------------------------------------------------
        uint s_i;
        if ( !parse_gph_op_lookahead( s_i, gph, gph_end ) ) break;

        //-------------------------------------------------------------------
        // Get op precedence and check that it's high enough.
        //-------------------------------------------------------------------
        uint op_prec = 0;
        if ( strcmp( &strings[s_i], "+" ) == 0 ||
             strcmp( &strings[s_i], "-" ) == 0 ) {
            op_prec = 10;
        } else if ( strcmp( &strings[s_i], "*" ) == 0 ||
                    strcmp( &strings[s_i], "/" ) == 0 ) {
            op_prec = 20;
        } else {
            die( "bad operator in expression: " + std::string( &strings[s_i] ) );
        }
        mdout << "op_prec=" << op_prec << " curr_prec=" << curr_prec << "\n";
        if ( op_prec < curr_prec ) break;
        curr_prec = op_prec;

        //-------------------------------------------------------------------
        // Consume operator and parse term.
        //-------------------------------------------------------------------
        uint op_i;
        if ( !parse_gph_op( op_i, gph, gph_end ) ) return false;

        uint opnd1_i;
        if ( !parse_gph_term( opnd1_i, gph, gph_end ) ) return false;

        //-------------------------------------------------------------------
        // Make an n-ary node out of this.
        // expr_i still contains opnd0_i
        //-------------------------------------------------------------------
        uint nary_i = gph_node_alloc( GRAPH_NODE_KIND::NARY );
        graph_nodes[nary_i].u.child_first_i = op_i;
        graph_nodes[op_i].parent_i = nary_i;
        graph_nodes[op_i].sibling_i = expr_i;
        graph_nodes[expr_i].parent_i = nary_i;
        graph_nodes[expr_i].sibling_i = opnd1_i;
        graph_nodes[opnd1_i].parent_i = nary_i;
        expr_i = nary_i;
    }

    return true;
}

inline bool Model::parse_real3( Model::real3& r3, const char *& xxx, const char * xxx_end, bool has_brackets )
{
    return (!has_brackets || expect_char( '[', xxx, xxx_end, true )) &&
           parse_real( r3.c[0], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[1], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[2], xxx, xxx_end, has_brackets ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Model::parse_real2( Model::real2& r2, const char *& xxx, const char * xxx_end, bool has_brackets )
{
    return (!has_brackets || expect_char( '[', xxx, xxx_end, true )) &&
           parse_real( r2.c[0], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r2.c[1], xxx, xxx_end, has_brackets ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Model::parse_real( Model::real& r, const char *& xxx, const char * xxx_end, bool skip_whitespace_first ) 
{
    Model::real64 r64;
    if ( !parse_real64( r64, xxx, xxx_end, skip_whitespace_first ) ) return false;
    r = r64;
    return true;
}

inline bool Model::parse_real64( Model::real64& r64, const char *& xxx, const char * xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );   // can span lines unlike below
    std::string s = "";
    bool in_frac = false;
    bool has_exp = false;
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    while( xxx != xxx_end )
    {
        char ch = *xxx;

        if ( ch == 'n' || ch == 'N' ) {
            // better be a NaN
            xxx++;
            if ( xxx == xxx_end || (*xxx != 'a' && *xxx != 'A') ) return false;
            xxx++;
            if ( xxx == xxx_end || (*xxx != 'n' && *xxx != 'N') ) return false;
            xxx++;
            r64 = 0.0;                    // make them zeros
            return true;
        }

        if ( ch == '-' && !in_frac ) {
            s += "-";
            xxx++;
            continue;
        }

        if ( ch == '.' && !in_frac ) {
            s += ".";
            in_frac = true;
            xxx++;
            continue;
        }

        if ( ch == 'e' || ch == 'E' ) {
            rtn_assert( !has_exp, "real has more than one 'e' exponent" );
            has_exp = true;
            s += std::string( 1, ch );
            xxx++;
            _int e10;
            if ( !parse_int( e10, xxx, xxx_end ) ) return false;
            s += std::to_string( e10 );
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;

        s += std::string( 1, ch );
        xxx++;
    }

    rtn_assert( s.length() != 0, "unable to parse real in file " + surrounding_lines( xxx, xxx_end ) );

    r64 = std::atof( s.c_str() );
    dprint( "real=" + std::to_string( r ) );
    return true;
}

inline bool Model::parse_int64( _int64& i, const char *& xxx, const char * xxx_end )
{
    bool vld = false;
    i = 0;
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    bool is_neg = false;
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( ch == '-' ) {
            rtn_assert( !is_neg, "too many minus signs" );
            is_neg = true;
            xxx++;
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;
        xxx++;

        i = i*10 + (ch - '0');
        vld = true;
    }

    if ( is_neg ) i = -i;
    rtn_assert( vld, "unable to parse int" + surrounding_lines( xxx, xxx_end ) );
    return true;
}

inline bool Model::parse_int( _int& i, const char *& xxx, const char * xxx_end )
{
    _int64 i64;
    bool r = parse_int64( i64, xxx, xxx_end );
    i = i64;
    return r;
}

inline bool Model::parse_uint64( uint64& u, const char *& xxx, const char * xxx_end, uint base )
{
    bool vld = false;
    u = 0;
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    while( xxx != xxx_end )
    {
        char ch = *xxx;

        if ( base == 10 ) {
            if ( ch >= '0' && ch <= '9' ) {
                u = u*10 + (ch - '0');
            } else {
                break;
            }
        } else if ( base == 16 ) {
            if ( ch >= '0' && ch <= '9' ) {
                u = u*16 + (ch - '0');
            } else if ( ch >= 'a' && ch <= 'f' ) {
                u = u*16 + (ch - 'a' + 10);
            } else if ( ch >= 'A' && ch <= 'F' ) {
                u = u*16 + (ch - 'A' + 10);
            } else {
                break;
            }
        } else if ( base == 8 ) {
            if ( ch >= '0' && ch <= '7' ) {
                u = u*8 + (ch - '0');
            } else {
                break;
            }
        } else { 
            rtn_assert( false, "bad base, must be 8, 10, or 16; 0 not yet supported" );
        }
        xxx++;
        vld = true;
    }

    rtn_assert( vld, "unable to parse uint64" + surrounding_lines( xxx, xxx_end ) );
    return true;
}

inline bool Model::parse_uint( uint& u, const char *& xxx, const char * xxx_end, uint base )
{
    uint64 u64;
    bool r = parse_uint64( u64, xxx, xxx_end, base );
    u = u64;
    return r;
}

inline bool Model::parse_bool( bool& b, const char *& xxx, const char * xxx_end )
{
    std::string id;
    if ( !parse_id( id , xxx, xxx_end ) ) return false;
    b = id == std::string( "true" );
    return true;
}

std::string Model::surrounding_lines( const char *& xxx, const char * xxx_end )
{
    uint eol_cnt = 0;
    std::string s = "";
    while( eol_cnt != 10 && xxx != xxx_end )
    {
        s += std::string( 1, *xxx ) ;
        if ( *xxx == '\n' ) eol_cnt++;
        xxx++;
    }
    return s;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//
// ASTC TEXTURE COMPRESSION AND DECOMPRESSION
//
// ASTC is a modern texture format that achieves excellent compression while
// maintaining image quality.
//
// ARM has an open-source compression and decompression program that acts as
// the golden standard implementation.
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
inline Model::real4 Model::Texture::texel_read_astc( const Model * model, uint mip_level, uint64 ui, uint64 vi, uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const
{
    //---------------------------------------------------------------
    // Each mipmap level has its own ASTC_Header and thus potential block footprint.
    // Find the beginning of the appropriate mip texture data.
    // All blocks are 16B, but texels represented can differ.
    //---------------------------------------------------------------
    const ASTC_Header * astc_hdr = reinterpret_cast<const ASTC_Header *>( &model->texels[texel_i] );
    die_assert( astc_hdr->magic_is_good(), "ASTC mip 0 header is corrupted - bad magic number" );
    uint mw = width;
    uint mh = height;
    if ( model->hdr->mipmap_filter != MIPMAP_FILTER::NONE ) {
        for( _int imip_level = mip_level; imip_level > 0 && !(mw == 1 && mh == 1); imip_level-- ) 
        {
            die_assert( astc_hdr->magic_is_good(), "ASTC mip " + std::to_string(imip_level) + " header is corrupted - bad magic number" );
            die_assert( astc_hdr->size_x() == mw,  "ASTC mip " + std::to_string(imip_level) + " header is corrupted - bad width" );
            die_assert( astc_hdr->size_y() == mh,  "ASTC mip " + std::to_string(imip_level) + " header is corrupted - bad height" );
            die_assert( astc_hdr->size_z() == 1,   "ASTC mip " + std::to_string(imip_level) + " header is corrupted - depth != 1" );

            if ( mw != 1 ) mw >>= 1;
            if ( mh != 1 ) mh >>= 1;
            astc_hdr += astc_hdr->blk_cnt();
        }
    }
    const unsigned char * bdata = reinterpret_cast<const unsigned char *>( astc_hdr+1 );

    //---------------------------------------------------------------
    // Find block x,y and offsets within block.
    // Find block data.
    //---------------------------------------------------------------
    uint bx = ui / astc_hdr->blockdim_x;
    uint by = vi / astc_hdr->blockdim_y;
    uint s  = ui - bx*astc_hdr->blockdim_x;
    uint t  = vi - by*astc_hdr->blockdim_y;
    uint r  = 0;
    mdout << "\nastc: uv=[" << ui << "," << vi << "] bxy=[" << bx << "," << by << "] st=[" << s << "," << t << "]\n";
    uint bw = astc_hdr->blk_cnt_x();
    bdata += 16 * (bx + by*bw);

    //---------------------------------------------------------------
    // Decode ASTC texel at [s,t,r] within the block.
    //---------------------------------------------------------------
    real4 rgba = astc_decode_texel( bdata, astc_hdr, s, t, r, do_srgb_to_linear );

    //---------------------------------------------------------------
    // Return texel and other info.
    //---------------------------------------------------------------
    if ( vaddr != nullptr )     *vaddr = reinterpret_cast<uint64_t>( bdata );
    if ( byte_cnt != nullptr )  *byte_cnt = 16;
    return rgba;
}

const Model::Texture::ASTC_Range_Encoding Model::Texture::astc_range_encodings[21] = 
{
    // these are used for both weights and colors
    {  2, 0, 0, 1 },    
    {  3, 1, 0, 0 },
    {  4, 0, 0, 2 },
    {  5, 0, 1, 0 },
    {  6, 1, 0, 1 },
    {  8, 0, 0, 3 },
    { 10, 0, 1, 1 },    
    { 12, 1, 0, 2 },
    { 16, 0, 0, 4 },
    { 20, 0, 1, 2 },
    { 24, 1, 0, 3 },
    { 32, 0, 0, 5 },
    { 40, 0, 1, 3 },
    { 48, 1, 0, 4 },
    { 64, 0, 0, 6 },
    { 80, 0, 1, 4 },
    { 96, 1, 0, 5 },
    {128, 0, 0, 7 },
    {160, 0, 1, 5 },
    {192, 1, 0, 6 },
    {256, 0, 0, 8 },
};

Model::Texture::ASTC_Decimation_Encoding Model::Texture::astc_weight_decimation_encodings;

const Model::uint Model::Texture::astc_color_quantization_mode[17][128] = 
{
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   // 255 == invalid
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    },
    {
        255, 255, 0, 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 0, 0, 0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 
        8, 8, 9, 10, 11, 11, 12, 13, 14, 14, 15, 16, 17, 17, 18, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 
        4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 
        12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 1, 1, 1, 
        2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 
        8, 8, 8, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 13, 13, 13, 
        14, 14, 14, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 19, 19, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 
        1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 
        5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 
        10, 10, 11, 11, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 
        15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 19, 19, 19, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 
        0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 
        4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 
        8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 
        12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 
        16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
        2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 
        6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 
        9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 
        13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 
        16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 
        5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 
        8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 
        11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 
        14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 
        17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 
        1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 
        4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 
        6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 
        9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 
        12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 
        14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 17, 17, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 
        3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 
        5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 
        8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 
        10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 
        13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
        2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 
        4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 
        7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 
        9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 
        11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 
        4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
        8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 
        10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 
        3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 
        5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 
        7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 
        8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 
        2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 
        4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 
        6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 
        5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 
    },
    {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 
        6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
    },
};
const Model::uint Model::Texture::astc_integer_trits[256][5] = 
{
    {0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}, {0, 0, 2, 0, 0},
    {0, 1, 0, 0, 0}, {1, 1, 0, 0, 0}, {2, 1, 0, 0, 0}, {1, 0, 2, 0, 0},
    {0, 2, 0, 0, 0}, {1, 2, 0, 0, 0}, {2, 2, 0, 0, 0}, {2, 0, 2, 0, 0},
    {0, 2, 2, 0, 0}, {1, 2, 2, 0, 0}, {2, 2, 2, 0, 0}, {2, 0, 2, 0, 0},
    {0, 0, 1, 0, 0}, {1, 0, 1, 0, 0}, {2, 0, 1, 0, 0}, {0, 1, 2, 0, 0},
    {0, 1, 1, 0, 0}, {1, 1, 1, 0, 0}, {2, 1, 1, 0, 0}, {1, 1, 2, 0, 0},
    {0, 2, 1, 0, 0}, {1, 2, 1, 0, 0}, {2, 2, 1, 0, 0}, {2, 1, 2, 0, 0},
    {0, 0, 0, 2, 2}, {1, 0, 0, 2, 2}, {2, 0, 0, 2, 2}, {0, 0, 2, 2, 2},
    {0, 0, 0, 1, 0}, {1, 0, 0, 1, 0}, {2, 0, 0, 1, 0}, {0, 0, 2, 1, 0},
    {0, 1, 0, 1, 0}, {1, 1, 0, 1, 0}, {2, 1, 0, 1, 0}, {1, 0, 2, 1, 0},
    {0, 2, 0, 1, 0}, {1, 2, 0, 1, 0}, {2, 2, 0, 1, 0}, {2, 0, 2, 1, 0},
    {0, 2, 2, 1, 0}, {1, 2, 2, 1, 0}, {2, 2, 2, 1, 0}, {2, 0, 2, 1, 0},
    {0, 0, 1, 1, 0}, {1, 0, 1, 1, 0}, {2, 0, 1, 1, 0}, {0, 1, 2, 1, 0},
    {0, 1, 1, 1, 0}, {1, 1, 1, 1, 0}, {2, 1, 1, 1, 0}, {1, 1, 2, 1, 0},
    {0, 2, 1, 1, 0}, {1, 2, 1, 1, 0}, {2, 2, 1, 1, 0}, {2, 1, 2, 1, 0},
    {0, 1, 0, 2, 2}, {1, 1, 0, 2, 2}, {2, 1, 0, 2, 2}, {1, 0, 2, 2, 2},
    {0, 0, 0, 2, 0}, {1, 0, 0, 2, 0}, {2, 0, 0, 2, 0}, {0, 0, 2, 2, 0},
    {0, 1, 0, 2, 0}, {1, 1, 0, 2, 0}, {2, 1, 0, 2, 0}, {1, 0, 2, 2, 0},
    {0, 2, 0, 2, 0}, {1, 2, 0, 2, 0}, {2, 2, 0, 2, 0}, {2, 0, 2, 2, 0},
    {0, 2, 2, 2, 0}, {1, 2, 2, 2, 0}, {2, 2, 2, 2, 0}, {2, 0, 2, 2, 0},
    {0, 0, 1, 2, 0}, {1, 0, 1, 2, 0}, {2, 0, 1, 2, 0}, {0, 1, 2, 2, 0},
    {0, 1, 1, 2, 0}, {1, 1, 1, 2, 0}, {2, 1, 1, 2, 0}, {1, 1, 2, 2, 0},
    {0, 2, 1, 2, 0}, {1, 2, 1, 2, 0}, {2, 2, 1, 2, 0}, {2, 1, 2, 2, 0},
    {0, 2, 0, 2, 2}, {1, 2, 0, 2, 2}, {2, 2, 0, 2, 2}, {2, 0, 2, 2, 2},
    {0, 0, 0, 0, 2}, {1, 0, 0, 0, 2}, {2, 0, 0, 0, 2}, {0, 0, 2, 0, 2},
    {0, 1, 0, 0, 2}, {1, 1, 0, 0, 2}, {2, 1, 0, 0, 2}, {1, 0, 2, 0, 2},
    {0, 2, 0, 0, 2}, {1, 2, 0, 0, 2}, {2, 2, 0, 0, 2}, {2, 0, 2, 0, 2},
    {0, 2, 2, 0, 2}, {1, 2, 2, 0, 2}, {2, 2, 2, 0, 2}, {2, 0, 2, 0, 2},
    {0, 0, 1, 0, 2}, {1, 0, 1, 0, 2}, {2, 0, 1, 0, 2}, {0, 1, 2, 0, 2},
    {0, 1, 1, 0, 2}, {1, 1, 1, 0, 2}, {2, 1, 1, 0, 2}, {1, 1, 2, 0, 2},
    {0, 2, 1, 0, 2}, {1, 2, 1, 0, 2}, {2, 2, 1, 0, 2}, {2, 1, 2, 0, 2},
    {0, 2, 2, 2, 2}, {1, 2, 2, 2, 2}, {2, 2, 2, 2, 2}, {2, 0, 2, 2, 2},
    {0, 0, 0, 0, 1}, {1, 0, 0, 0, 1}, {2, 0, 0, 0, 1}, {0, 0, 2, 0, 1},
    {0, 1, 0, 0, 1}, {1, 1, 0, 0, 1}, {2, 1, 0, 0, 1}, {1, 0, 2, 0, 1},
    {0, 2, 0, 0, 1}, {1, 2, 0, 0, 1}, {2, 2, 0, 0, 1}, {2, 0, 2, 0, 1},
    {0, 2, 2, 0, 1}, {1, 2, 2, 0, 1}, {2, 2, 2, 0, 1}, {2, 0, 2, 0, 1},
    {0, 0, 1, 0, 1}, {1, 0, 1, 0, 1}, {2, 0, 1, 0, 1}, {0, 1, 2, 0, 1},
    {0, 1, 1, 0, 1}, {1, 1, 1, 0, 1}, {2, 1, 1, 0, 1}, {1, 1, 2, 0, 1},
    {0, 2, 1, 0, 1}, {1, 2, 1, 0, 1}, {2, 2, 1, 0, 1}, {2, 1, 2, 0, 1},
    {0, 0, 1, 2, 2}, {1, 0, 1, 2, 2}, {2, 0, 1, 2, 2}, {0, 1, 2, 2, 2},
    {0, 0, 0, 1, 1}, {1, 0, 0, 1, 1}, {2, 0, 0, 1, 1}, {0, 0, 2, 1, 1},
    {0, 1, 0, 1, 1}, {1, 1, 0, 1, 1}, {2, 1, 0, 1, 1}, {1, 0, 2, 1, 1},
    {0, 2, 0, 1, 1}, {1, 2, 0, 1, 1}, {2, 2, 0, 1, 1}, {2, 0, 2, 1, 1},
    {0, 2, 2, 1, 1}, {1, 2, 2, 1, 1}, {2, 2, 2, 1, 1}, {2, 0, 2, 1, 1},
    {0, 0, 1, 1, 1}, {1, 0, 1, 1, 1}, {2, 0, 1, 1, 1}, {0, 1, 2, 1, 1},
    {0, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {2, 1, 1, 1, 1}, {1, 1, 2, 1, 1},
    {0, 2, 1, 1, 1}, {1, 2, 1, 1, 1}, {2, 2, 1, 1, 1}, {2, 1, 2, 1, 1},
    {0, 1, 1, 2, 2}, {1, 1, 1, 2, 2}, {2, 1, 1, 2, 2}, {1, 1, 2, 2, 2},
    {0, 0, 0, 2, 1}, {1, 0, 0, 2, 1}, {2, 0, 0, 2, 1}, {0, 0, 2, 2, 1},
    {0, 1, 0, 2, 1}, {1, 1, 0, 2, 1}, {2, 1, 0, 2, 1}, {1, 0, 2, 2, 1},
    {0, 2, 0, 2, 1}, {1, 2, 0, 2, 1}, {2, 2, 0, 2, 1}, {2, 0, 2, 2, 1},
    {0, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {2, 2, 2, 2, 1}, {2, 0, 2, 2, 1},
    {0, 0, 1, 2, 1}, {1, 0, 1, 2, 1}, {2, 0, 1, 2, 1}, {0, 1, 2, 2, 1},
    {0, 1, 1, 2, 1}, {1, 1, 1, 2, 1}, {2, 1, 1, 2, 1}, {1, 1, 2, 2, 1},
    {0, 2, 1, 2, 1}, {1, 2, 1, 2, 1}, {2, 2, 1, 2, 1}, {2, 1, 2, 2, 1},
    {0, 2, 1, 2, 2}, {1, 2, 1, 2, 2}, {2, 2, 1, 2, 2}, {2, 1, 2, 2, 2},
    {0, 0, 0, 1, 2}, {1, 0, 0, 1, 2}, {2, 0, 0, 1, 2}, {0, 0, 2, 1, 2},
    {0, 1, 0, 1, 2}, {1, 1, 0, 1, 2}, {2, 1, 0, 1, 2}, {1, 0, 2, 1, 2},
    {0, 2, 0, 1, 2}, {1, 2, 0, 1, 2}, {2, 2, 0, 1, 2}, {2, 0, 2, 1, 2},
    {0, 2, 2, 1, 2}, {1, 2, 2, 1, 2}, {2, 2, 2, 1, 2}, {2, 0, 2, 1, 2},
    {0, 0, 1, 1, 2}, {1, 0, 1, 1, 2}, {2, 0, 1, 1, 2}, {0, 1, 2, 1, 2},
    {0, 1, 1, 1, 2}, {1, 1, 1, 1, 2}, {2, 1, 1, 1, 2}, {1, 1, 2, 1, 2},
    {0, 2, 1, 1, 2}, {1, 2, 1, 1, 2}, {2, 2, 1, 1, 2}, {2, 1, 2, 1, 2},
    {0, 2, 2, 2, 2}, {1, 2, 2, 2, 2}, {2, 2, 2, 2, 2}, {2, 1, 2, 2, 2}
};

const Model::uint Model::Texture::astc_integer_quints[128][3] = 
{
    {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0},
    {4, 0, 0}, {0, 4, 0}, {4, 4, 0}, {4, 4, 4},
    {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {3, 1, 0},
    {4, 1, 0}, {1, 4, 0}, {4, 4, 1}, {4, 4, 4},
    {0, 2, 0}, {1, 2, 0}, {2, 2, 0}, {3, 2, 0},
    {4, 2, 0}, {2, 4, 0}, {4, 4, 2}, {4, 4, 4},
    {0, 3, 0}, {1, 3, 0}, {2, 3, 0}, {3, 3, 0},
    {4, 3, 0}, {3, 4, 0}, {4, 4, 3}, {4, 4, 4},
    {0, 0, 1}, {1, 0, 1}, {2, 0, 1}, {3, 0, 1},
    {4, 0, 1}, {0, 4, 1}, {4, 0, 4}, {0, 4, 4},
    {0, 1, 1}, {1, 1, 1}, {2, 1, 1}, {3, 1, 1},
    {4, 1, 1}, {1, 4, 1}, {4, 1, 4}, {1, 4, 4},
    {0, 2, 1}, {1, 2, 1}, {2, 2, 1}, {3, 2, 1},
    {4, 2, 1}, {2, 4, 1}, {4, 2, 4}, {2, 4, 4},
    {0, 3, 1}, {1, 3, 1}, {2, 3, 1}, {3, 3, 1},
    {4, 3, 1}, {3, 4, 1}, {4, 3, 4}, {3, 4, 4},
    {0, 0, 2}, {1, 0, 2}, {2, 0, 2}, {3, 0, 2},
    {4, 0, 2}, {0, 4, 2}, {2, 0, 4}, {3, 0, 4},
    {0, 1, 2}, {1, 1, 2}, {2, 1, 2}, {3, 1, 2},
    {4, 1, 2}, {1, 4, 2}, {2, 1, 4}, {3, 1, 4},
    {0, 2, 2}, {1, 2, 2}, {2, 2, 2}, {3, 2, 2},
    {4, 2, 2}, {2, 4, 2}, {2, 2, 4}, {3, 2, 4},
    {0, 3, 2}, {1, 3, 2}, {2, 3, 2}, {3, 3, 2},
    {4, 3, 2}, {3, 4, 2}, {2, 3, 4}, {3, 3, 4},
    {0, 0, 3}, {1, 0, 3}, {2, 0, 3}, {3, 0, 3},
    {4, 0, 3}, {0, 4, 3}, {0, 0, 4}, {1, 0, 4},
    {0, 1, 3}, {1, 1, 3}, {2, 1, 3}, {3, 1, 3},
    {4, 1, 3}, {1, 4, 3}, {0, 1, 4}, {1, 1, 4},
    {0, 2, 3}, {1, 2, 3}, {2, 2, 3}, {3, 2, 3},
    {4, 2, 3}, {2, 4, 3}, {0, 2, 4}, {1, 2, 4},
    {0, 3, 3}, {1, 3, 3}, {2, 3, 3}, {3, 3, 3},
    {4, 3, 3}, {3, 4, 3}, {0, 3, 4}, {1, 3, 4}
};

const Model::uint Model::Texture::astc_weight_unquantized[12][32] = 
{
    {0, 64},
    {0, 32, 64},
    {0, 21, 43, 64},
    {0, 16, 32, 48, 64},
    {0, 64, 12, 52, 25, 39},
    {0, 9, 18, 27, 37, 46, 55, 64},
    {0, 64, 7, 57, 14, 50, 21, 43, 28, 36},
    {0, 64, 17, 47, 5, 59, 23, 41, 11, 53, 28, 36},
    {0, 4, 8, 12, 17, 21, 25, 29, 35, 39, 43, 47, 52, 56, 60, 64},
    {0, 64, 16, 48, 3, 61, 19, 45, 6, 58, 23, 41, 9, 55, 26, 38, 13, 51, 29, 35},
    {0, 64, 8, 56, 16, 48, 24, 40, 2, 62, 11, 53, 19, 45, 27, 37, 5, 59, 13, 51, 22, 42, 30, 34},
    {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64},
};

const Model::uint Model::Texture::astc_color_unquantized[21][256] = 
{
    {
        0, 255
    },
    {
        0, 128, 255
    },
    {
        0, 85, 170, 255
    },
    {
        0, 64, 128, 192, 255
    },
    {
        0, 255, 51, 204, 102, 153
    },
    {
        0, 36, 73, 109, 146, 182, 219, 255
    },
    {
        0, 255, 28, 227, 56, 199, 84, 171, 113, 142
    },
    {
        0, 255, 69, 186, 23, 232, 92, 163, 46, 209, 116, 139
    },
    {
        0, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255
    },
    {
        0, 255, 67, 188, 13, 242, 80, 175, 27, 228, 94, 161, 40, 215, 107, 148,
        54, 201, 121, 134
    },
    {
        0, 255, 33, 222, 66, 189, 99, 156, 11, 244, 44, 211, 77, 178, 110, 145,
        22, 233, 55, 200, 88, 167, 121, 134
    },
    {
        0, 8, 16, 24, 33, 41, 49, 57, 66, 74, 82, 90, 99, 107, 115, 123,
        132, 140, 148, 156, 165, 173, 181, 189, 198, 206, 214, 222, 231, 239, 247, 255
    },
    {
        0, 255, 32, 223, 65, 190, 97, 158, 6, 249, 39, 216, 71, 184, 104, 151,
        13, 242, 45, 210, 78, 177, 110, 145, 19, 236, 52, 203, 84, 171, 117, 138,
        26, 229, 58, 197, 91, 164, 123, 132
    },
    {
        0, 255, 16, 239, 32, 223, 48, 207, 65, 190, 81, 174, 97, 158, 113, 142,
        5, 250, 21, 234, 38, 217, 54, 201, 70, 185, 86, 169, 103, 152, 119, 136,
        11, 244, 27, 228, 43, 212, 59, 196, 76, 179, 92, 163, 108, 147, 124, 131
    },
    {
        0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
        65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113, 117, 121, 125,
        130, 134, 138, 142, 146, 150, 154, 158, 162, 166, 170, 174, 178, 182, 186, 190,
        195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243, 247, 251, 255
    },
    {
        0, 255, 16, 239, 32, 223, 48, 207, 64, 191, 80, 175, 96, 159, 112, 143,
        3, 252, 19, 236, 35, 220, 51, 204, 67, 188, 83, 172, 100, 155, 116, 139,
        6, 249, 22, 233, 38, 217, 54, 201, 71, 184, 87, 168, 103, 152, 119, 136,
        9, 246, 25, 230, 42, 213, 58, 197, 74, 181, 90, 165, 106, 149, 122, 133,
        13, 242, 29, 226, 45, 210, 61, 194, 77, 178, 93, 162, 109, 146, 125, 130
    },
    {
        0, 255, 8, 247, 16, 239, 24, 231, 32, 223, 40, 215, 48, 207, 56, 199,
        64, 191, 72, 183, 80, 175, 88, 167, 96, 159, 104, 151, 112, 143, 120, 135,
        2, 253, 10, 245, 18, 237, 26, 229, 35, 220, 43, 212, 51, 204, 59, 196,
        67, 188, 75, 180, 83, 172, 91, 164, 99, 156, 107, 148, 115, 140, 123, 132,
        5, 250, 13, 242, 21, 234, 29, 226, 37, 218, 45, 210, 53, 202, 61, 194,
        70, 185, 78, 177, 86, 169, 94, 161, 102, 153, 110, 145, 118, 137, 126, 129
    },
    {
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
        32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
        64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94,
        96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126,
        129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159,
        161, 163, 165, 167, 169, 171, 173, 175, 177, 179, 181, 183, 185, 187, 189, 191,
        193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223,
        225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247, 249, 251, 253, 255
    },
    {
        0, 255, 8, 247, 16, 239, 24, 231, 32, 223, 40, 215, 48, 207, 56, 199,
        64, 191, 72, 183, 80, 175, 88, 167, 96, 159, 104, 151, 112, 143, 120, 135,
        1, 254, 9, 246, 17, 238, 25, 230, 33, 222, 41, 214, 49, 206, 57, 198,
        65, 190, 73, 182, 81, 174, 89, 166, 97, 158, 105, 150, 113, 142, 121, 134,
        3, 252, 11, 244, 19, 236, 27, 228, 35, 220, 43, 212, 51, 204, 59, 196,
        67, 188, 75, 180, 83, 172, 91, 164, 99, 156, 107, 148, 115, 140, 123, 132,
        4, 251, 12, 243, 20, 235, 28, 227, 36, 219, 44, 211, 52, 203, 60, 195,
        68, 187, 76, 179, 84, 171, 92, 163, 100, 155, 108, 147, 116, 139, 124, 131,
        6, 249, 14, 241, 22, 233, 30, 225, 38, 217, 46, 209, 54, 201, 62, 193,
        70, 185, 78, 177, 86, 169, 94, 161, 102, 153, 110, 145, 118, 137, 126, 129
    },
    {
        0, 255, 4, 251, 8, 247, 12, 243, 16, 239, 20, 235, 24, 231, 28, 227,
        32, 223, 36, 219, 40, 215, 44, 211, 48, 207, 52, 203, 56, 199, 60, 195,
        64, 191, 68, 187, 72, 183, 76, 179, 80, 175, 84, 171, 88, 167, 92, 163,
        96, 159, 100, 155, 104, 151, 108, 147, 112, 143, 116, 139, 120, 135, 124, 131,
        1, 254, 5, 250, 9, 246, 13, 242, 17, 238, 21, 234, 25, 230, 29, 226,
        33, 222, 37, 218, 41, 214, 45, 210, 49, 206, 53, 202, 57, 198, 61, 194,
        65, 190, 69, 186, 73, 182, 77, 178, 81, 174, 85, 170, 89, 166, 93, 162,
        97, 158, 101, 154, 105, 150, 109, 146, 113, 142, 117, 138, 121, 134, 125, 130,
        2, 253, 6, 249, 10, 245, 14, 241, 18, 237, 22, 233, 26, 229, 30, 225,
        34, 221, 38, 217, 42, 213, 46, 209, 50, 205, 54, 201, 58, 197, 62, 193,
        66, 189, 70, 185, 74, 181, 78, 177, 82, 173, 86, 169, 90, 165, 94, 161,
        98, 157, 102, 153, 106, 149, 110, 145, 114, 141, 118, 137, 122, 133, 126, 129
    },
    {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
        48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
        64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
        80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
        96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
        112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
    }
};

static inline uint64_t _astc_bits( const unsigned char * bdata, uint start, uint len ) 
{
    // extract bits from 128 bits
    uint64_t bits_lo = *reinterpret_cast<const uint64_t *>( bdata+0 );
    uint64_t bits_hi = *reinterpret_cast<const uint64_t *>( bdata+8 );
    uint64_t bits = bits_lo >> start;
    if ( (start+len) <= 64 ) {
        // lo only
        bits &= (1ULL << len) - 1ULL;
    } else if ( start < 64 ) {
        // mixed
        uint64_t lower_len = 64 - start;
        uint64_t upper_len = len - lower_len;
        uint64_t upper     = bits_hi & ((1ULL << upper_len)-1);
        bits |= upper << lower_len;
    } else {
        // hi only
        bits = (bits_hi >> (start-64)) & ((1ULL << len)-1);
    }
    return bits;
}

inline void Model::Texture::astc_blk_dims_get( uint width, uint height, uint& blk_dim_x, uint& blk_dim_y )
{
    //---------------------------------------------------------------
    // Always choose squares for now using max(width, height).
    // Algorithm here is arbitrary.
    //---------------------------------------------------------------
    uint dim = std::max( width, height );
         if ( dim >= 2048 )     blk_dim_x = 12;
    else if ( dim >=  512 )     blk_dim_x = 10;
    else if ( dim >=  256 )     blk_dim_x =  8;
    else if ( dim >=   64 )     blk_dim_x =  6;
    else if ( dim >=   32 )     blk_dim_x =  5;
    else                        blk_dim_x =  4;
    blk_dim_y = blk_dim_x;
}

void Model::Texture::astc_blk_addr_get( uint      ray_id,  
                                        uint      blockdim_x, uint blockdim_y, uint blockdim_z,
                                        uint      xsize,      uint ysize,      uint zsize,
                                        real      u,          real v,          real w,
                                        real      sqrt_frac_uv_covg, 
                                        uint      size_w,
                                        uint      addr_w,
                                        uint      uv_frac_w,
                                        uint      covg_frac_w,
                                        uint64_t  blks_addr,
                                        uint64_t& blk_addr, uint& s, uint& t, uint& r )
{
    const uint64_t sqrt_xsize_ysize_zsize = std::sqrt( real(xsize) * real(ysize) * real(zsize) );
    const uint64_t u_fxd = u * real( 1 << uv_frac_w );   
    const uint64_t v_fxd = v * real( 1 << uv_frac_w );
    const uint64_t w_fxd = w * real( 1 << uv_frac_w );
    const uint64_t sqrt_frac_uv_covg_fxd = int( sqrt_frac_uv_covg * real(1 << covg_frac_w) );

    mdout << "astc_blk_addr_get: ray_id=" << ray_id <<
                               " blockdim_x=" << blockdim_x << " blockdim_y=" << blockdim_y << " blockdim_z=" << blockdim_z <<
                               std::hex << " xsize=0x" << xsize << " ysize=0x" << ysize << " zsize=0x" << zsize << 
                               " sqrt_xsize_ysize_zsize=0x" << sqrt_xsize_ysize_zsize << 
                               std::dec << " u=" << u << " v=" << v << " w=" << w << 
                               std::hex << " u_fxd=0x" << u_fxd << " v_fxd=0x" << v_fxd << " w_fxd=0x" << w_fxd << 
                               " sqrt_frac_uv_covg_fxd=0x" << sqrt_frac_uv_covg_fxd <<
                               std::dec << " size_w=" << size_w << " addr_w=" << addr_w << " uv_frac_w=" << uv_frac_w << 
                               " covg_frac_w=" << covg_frac_w << 
                               std::hex << " blks_addr=0x" << blks_addr << std::dec << "\n";

    //---------------------------------------------------------------
    // MIP LEVEL
    // Find the proper mip level.  We compute which mip levels are viable, in parallel,' )
    // then use a priority encoder to pixel the first viable one.' )
    // Mip level 0 is the finest level of the texture.' )
    //---------------------------------------------------------------
    const uint64_t iwidth_of_footprint = (sqrt_frac_uv_covg_fxd * sqrt_xsize_ysize_zsize) >> covg_frac_w; // integer width
    const uint64_t iwidth_of_footprint_lg2 = (iwidth_of_footprint == 0) ? 0 : uint( std::log2( real(iwidth_of_footprint) ) );
    uint64_t mip_xsize = xsize;
    uint64_t mip_ysize = ysize;
    uint64_t mip_zsize = zsize;
    uint64_t mip_blk_cnt_x = 0;
    uint64_t mip_blk_cnt_y = 0;
    uint64_t mip_blk_cnt_z = 0;
    uint64_t mip_blk_cnt = 0;
    uint64_t mip_blk0_offset = 0;
    uint i = 0;
    for( ; i < size_w; i++ )
    {
        mip_blk_cnt_x = (mip_xsize + blockdim_x - 1) / blockdim_x;
        mip_blk_cnt_y = (mip_ysize + blockdim_y - 1) / blockdim_y;
        mip_blk_cnt_z = (mip_zsize + blockdim_z - 1) / blockdim_z;
        mip_blk_cnt   = mip_blk_cnt_x * mip_blk_cnt_y * mip_blk_cnt_z;

        mdout << "astc_blk_addr_get: mip" << i << ": ray_id=" << ray_id << 
                std::hex << " mip_blk0_offset=0x" << mip_blk0_offset << 
                " mip_blk_cnt_x=0x" << mip_blk_cnt_x << " mip_blk_cnt_y=0x" << mip_blk_cnt_y << " mip_blk_cnt_z=0x" << mip_blk_cnt_z << 
                " mip_blk_cnt=0x" << mip_blk_cnt << 
                std::dec << "\n";

        if ( i == (size_w-1) || i >= iwidth_of_footprint_lg2 ) break;

        mip_blk0_offset += mip_blk_cnt + 1; // skip the next header, too
        mip_xsize >>= 1;
        mip_ysize >>= 1;
        mip_zsize >>= 1;
        if ( mip_xsize == 0 ) mip_xsize = 1;
        if ( mip_ysize == 0 ) mip_ysize = 1;
        if ( mip_zsize == 0 ) mip_zsize = 1;
    }

    mdout << "astc_blk_addr_get: ray_id=" << ray_id << 
            std::hex << " iwidth_of_footprint=0x" << iwidth_of_footprint << std::dec <<
            " iwidth_of_footprint_lg2=" << iwidth_of_footprint_lg2 << " mip_i=" << i <<
            std::hex << " mip_xsize=0x" << mip_xsize << " mip_ysize=0x" << mip_ysize << " mip_zsize=0x" << mip_zsize << 
            std::hex << " mip_blk_cnt_x=0x" << mip_blk_cnt_x << " mip_blk_cnt_y=0x" << mip_blk_cnt_y << " mip_blk_cnt_z=0x" << mip_blk_cnt_z << 
            std::hex << " mip_blk_cnt=0x" << mip_blk_cnt << " mip_blk0_offset=0x" << mip_blk0_offset <<
            std::dec << "\n";

    //---------------------------------------------------------------' )
    // Find block bx,by within map texture, block address, and texel offsets s,t within block.' )
    //---------------------------------------------------------------' )
    uint64_t ui = (u_fxd * mip_xsize) >> uv_frac_w;
             ui %= mip_xsize;
    uint64_t vi = (v_fxd * mip_ysize) >> uv_frac_w;
             vi %= mip_ysize;
    uint64_t wi = (w_fxd * mip_zsize) >> uv_frac_w;
             wi %= mip_zsize;
    uint64_t bx = ui / blockdim_x;
    uint64_t by = vi / blockdim_y;
    uint64_t bz = wi / blockdim_z;

    mdout << "astc_blk_addr_get: ray_id=" << ray_id << 
                        std::hex << " u_fxd=0x" << u_fxd << " v_fxd=0x" << v_fxd << " w_fxd=0x" << w_fxd <<
                        " ui=0x" << ui << " vi=0x" << vi << " wi=0x" << wi <<
                        " bx=0x" << bx << " by=0x" << by << " bz=0x" << bz << std::dec << "\n";

    blk_addr = blks_addr + mip_blk0_offset + bx + by*mip_blk_cnt_x + bz*mip_blk_cnt_x*mip_blk_cnt_y;
    s = ui - bx*blockdim_x;
    t = vi - by*blockdim_y;
    r = wi - bz*blockdim_z;

    mdout << "astc_blk_addr_get: ray_id=" << ray_id <<
                               std::hex << " blk_addr=0x" << blk_addr << std::dec <<
                               " s=" << s << " t=" << t << " r=" << r << "\n";

    die_assert( s < blockdim_x, "s >= blockdim_x" );
    die_assert( t < blockdim_y, "t >= blockdim_y" );
    die_assert( r < blockdim_z, "r >= blockdim_z" );
}

#define _bits(  start, len ) _astc_bits( bdata,  start, len )
#define _rbits( start, len ) _astc_bits( rbdata, start, len )

inline std::string Model::Texture::astc_dat_str( const unsigned char * bdata ) 
{
    static const char * digits[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f" };

    std::string s = "";
    for( uint i = 0; i < 16; i++ )
    {
        uint b = bdata[i];
        for( uint j = 0; j < 2; j++ )
        {
            uint n = (b >> (j*4)) & 0xf;
            s = digits[n] + s;
        }
    }
    return "0x" + s;
}

Model::real4 Model::Texture::astc_decode_texel( const unsigned char * bdata, const ASTC_Header * astc_hdr, uint s, uint t, uint r, bool do_srgb_to_linear )
{
    mdout << "astc: " << *astc_hdr << " s=" << s << " t=" << t << " r=" << r << 
             " do_srgb_to_linear=" << do_srgb_to_linear << " dat=" << astc_dat_str( bdata ) << "\n";

    //---------------------------------------------------------------
    // Determine if void-extent.
    //---------------------------------------------------------------
    real4 rgba;
    uint block_mode = _bits(0, 11);
    if ( (block_mode & 0x1ff) == 0x1fc ) {
        //---------------------------------------------------------------
        // VOID EXTENT CASE
        //
        // Format can be floating-point or unorm16.
        // Color bits are [127:64].
        // Return constant color.
        //---------------------------------------------------------------
        bool is_hdr = (block_mode & 0x200) != 0;
        die_assert( !is_hdr, "ASTC can't currently handle HDR void-extent blocks" );
        for( uint c = 0; c < 4; c++ )
        {
            bool do_srgb = do_srgb_to_linear && c != 3;
            uint C = _bits(64 + 16*c, 16);
            mdout << "C[" << c << "]=" << std::hex << C << std::dec << "\n";
            if ( do_srgb ) {
                C = (C >> 8) & 0xff;
                rgba.c[c] = srgb8_to_linear_gamma( C );
            } else {
                rgba.c[c] = (C == 0xffff) ? 1.0 : (real(C) / real(1 << 16));
            }
        }
        mdout << "astc: void extent rgba=" << rgba << "\n";
        return rgba;
    }

    //---------------------------------------------------------------
    // NORMAL CASE
    //
    // Refer to the ASTC Specification 1.0, which we follow here.
    // rbdata == reversal of 128 bdata bits so that we can easily traverse from other end.
    // plane_cnt == 1 + D bit from spec
    // partition_cnt
    // See Table 11 for block mode.
    //---------------------------------------------------------------
    uint texel_cnt = astc_hdr->blockdim_x * astc_hdr->blockdim_y * astc_hdr->blockdim_z;
    unsigned char rbdata[16];
    uint plane_cnt;
    uint partition_cnt;
    uint weights_w;  // M
    uint weights_h;  // N
    uint weights_d;  // Q
    uint R; // R = range encoding
    uint H; // H = high-precision?
    uint weights_qmode;

    astc_decode_block_mode( bdata, rbdata, plane_cnt, partition_cnt, weights_w, weights_h, weights_d, R, H, weights_qmode );

    uint plane_weights_cnt = weights_w * weights_h * weights_d;  // per-plane
    uint weights_cnt = plane_weights_cnt * plane_cnt;            // total

    mdout << "astc: blockdims=[" << uint32_t(astc_hdr->blockdim_x) << "," << uint32_t(astc_hdr->blockdim_y) << "," << uint32_t(astc_hdr->blockdim_z) << "]" <<
             " plane_cnt=" << plane_cnt << " partition_cnt=" << partition_cnt << 
             " weights_dims=[" << weights_w << "," << weights_h << "," << weights_d << "] R=" << R << " H=" << H <<
             " plane_weights_cnt=" << plane_weights_cnt << " weights_cnt=" << weights_cnt << "\n";

    //---------------------------------------------------------------
    // Use lookup table based on weights_qmode, which is derived from R and H.
    // Retrive trits, quints, and bits.
    // Color endpoints use the same table, which is a superset of our needs.
    //---------------------------------------------------------------
    die_assert( R >= 2, "astc: bad R range encoding" );
    const ASTC_Range_Encoding& enc = astc_range_encodings[weights_qmode];
    uint weight_trits  = enc.trits;
    uint weight_quints = enc.quints;
    uint weight_bits   = enc.bits;
    mdout << "astc: weight qmode=" << weights_qmode << " max_p1=" << enc.max_p1 << " trits=" << weight_trits << " quints=" << weight_quints << " bits=" << weight_bits << "\n";

    //---------------------------------------------------------------
    // See section 3.16.
    // Calculate some size info we need.
    //---------------------------------------------------------------
    uint weights_bit_cnt = (weights_cnt*8*weight_trits  + 4) / 5 +
                           (weights_cnt*7*weight_quints + 2) / 3 + 
                           (weights_cnt*1*weight_bits   + 0) / 1;

    uint weights_start       = 0;                       // counting from msb
    uint below_weights_start = 128 - weights_bit_cnt;   // counting from lsb

    mdout << "astc: weights_bit_cnt=" << weights_bit_cnt << " below_weights_start=" << below_weights_start << "\n";

    //---------------------------------------------------------------
    // Decode partition for our texel.
    //---------------------------------------------------------------
    uint partition = astc_decode_partition( bdata, s, t, r, partition_cnt, texel_cnt );

    //---------------------------------------------------------------
    // Determine the Color Endpoint Mode (CEM) for our partition.
    // See Table 13 and 14.
    // Also determine the number of trits, quints, and bits,
    // factoring in the space left.
    // See Table 15 and section 3.7.
    //---------------------------------------------------------------
    uint cem;
    uint cem_qmode;
    uint cem_trits;
    uint cem_quints;
    uint cem_bits;
    uint color_endpoints_cnt;
    uint color_endpoints_start;
    uint color_endpoints_v_first;  // for our partition
    uint color_endpoints_v_cnt;    // for our partition
    astc_decode_color_endpoint_mode( bdata, plane_cnt, partition, partition_cnt, weights_bit_cnt, below_weights_start, 
                                     cem, cem_qmode, cem_trits, cem_quints, cem_bits, 
                                     color_endpoints_cnt, color_endpoints_start, color_endpoints_v_first, color_endpoints_v_cnt );
    mdout << "astc: cem=" << cem << " cem_qmode=" << cem_qmode << " trits=" << cem_trits << " quints=" << cem_quints << " bits=" << cem_bits << 
            " color_endpoints_cnt=" << color_endpoints_cnt << " color_endpoints_start=" << color_endpoints_start << 
            " color_endpoints_v_first=" << color_endpoints_v_first << " color_endpoints_v_cnt=" << color_endpoints_v_cnt << "\n";

    //---------------------------------------------------------------
    // Decode the endpoints for our partition.
    //---------------------------------------------------------------
    uint color_endpoints[4][2];
    astc_decode_color_endpoints( bdata, cem, cem_qmode, cem_trits, cem_quints, cem_bits, 
                                 color_endpoints_cnt, color_endpoints_start, color_endpoints_v_first, color_endpoints_v_cnt,
                                 color_endpoints );

    //---------------------------------------------------------------
    // If dual-plane, find which channel is the dual-plane one.
    // There can be only one.
    //---------------------------------------------------------------
    uint plane1_chan = 0;
    if ( plane_cnt == 2 ) {
        uint plane1_chan_start = below_weights_start - 2;
        die_assert( plane1_chan_start < 127, "ASTC bad plane1_chan_start" );
        plane1_chan = _bits(plane1_chan_start, 2); 
    }

    //---------------------------------------------------------------
    // Get 1 or 2 plane weights.
    //---------------------------------------------------------------
    uint plane_weights[2];
    astc_decode_weights( rbdata, weights_start, plane_cnt, 
                         astc_hdr->blockdim_x, astc_hdr->blockdim_y, astc_hdr->blockdim_z,
                         weights_w, weights_h, weights_d,
                         s, t, r,
                         weight_trits, weight_quints, weight_bits, weights_qmode,
                         plane_weights );
    mdout << "astc: plane_weights[0]=" << plane_weights[0];
    if ( plane_cnt == 2 ) mdout << " plane_weights[1]=" << plane_weights[1];
    mdout << "\n";

    //---------------------------------------------------------------
    // See section 3.13.
    // Apply weights to color endpoints.
    // Assume non-sRGB for now.
    //---------------------------------------------------------------
    for( uint c = 0; c < 4; c++ )
    {
        uint pi = (plane_cnt == 2 && plane1_chan == c) ? 1 : 0;
        uint w  = plane_weights[pi];

        bool do_srgb = do_srgb_to_linear && c != 3;
        uint C0 = (color_endpoints[c][0] << 8) | (do_srgb ? 0x80 : color_endpoints[c][0]);
        uint C1 = (color_endpoints[c][1] << 8) | (do_srgb ? 0x80 : color_endpoints[c][1]);
        uint C = (C0*(64-w) + C1*w + 32) / 64;
        if ( do_srgb ) {
            C = (C >> 8) & 0xff;
            rgba.c[c] = srgb8_to_linear_gamma( C );
        } else {
            rgba.c[c] = (C == 0xffff) ? 1.0 : (real(C) / real(1 << 16));
        }
        uint rgba_fp = rgba.c[c] * 0x10000;
        mdout << "astc: c=" << c << " endpoint0=" << color_endpoints[c][0] << " endpoint1=" << color_endpoints[c][1] <<
                        std::hex << " C0=0x" << C0 << " C1=0x" << C1 << 
                        std::dec << " plane1_chan=" << plane1_chan << " pi=" << pi << " w=" << w << 
                        std::hex << " C=0x" << C << " rgba_fp[c]=0x" << rgba_fp << std::dec << " rgba[c]=" << rgba.c[c] << "\n";
    }
    return rgba;
}

inline void Model::Texture::astc_decode_block_mode( const unsigned char * bdata, unsigned char * rbdata, 
                                                    uint& plane_cnt, uint& partition_cnt, 
                                                    uint& weights_w, uint& weights_h, uint& weights_d, 
                                                    uint& R, uint& H, uint& weights_qmode )
{
    //---------------------------------------------------------------
    // Refer to the ASTC Specification 1.0, which we follow here.
    // First reverse bdata bits to make it easier to traverse from the other end.
    // plane_cnt == 1 + D bit from spec
    // partition_cnt
    // See Table 11 for block mode.
    //---------------------------------------------------------------
    for( uint i = 0; i < 16; i++ )
    {
        unsigned char b = bdata[15-i];
        b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;  // reverse bits in this byte
        b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
        b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
        rbdata[i] = b;
    }

    plane_cnt     = 1 + _bits(10, 1);    // hammered to 0 below for A6_B6 mode 
    partition_cnt = _bits(11, 2) + 1;

    uint a = _bits(5, 2);  
    uint b = _bits(7, 2); 
    R = _bits(4, 1);  // base qmode
    H = _bits(9, 1);  // high-precision?
    weights_d = 1;
    if ( _bits(0, 2) != 0 ) {
        // first 5 rows in table 11
        R |= _bits(0, 2) << 1;      
        uint mode_bits = _bits(2, 2);
        if ( mode_bits == 0 ) {
            // B4_A2
            weights_w = b + 4;
            weights_h = a + 2;
        
        } else if ( mode_bits == 1 ) {
            // B8_A2
            weights_w = b + 8;
            weights_h = a + 2;

        } else if ( mode_bits == 2 ) {
            // A2_B8
            weights_w = a + 2;
            weights_h = b + 8;

        } else if ( _bits(8, 1) ) {
            // B2_A2
            weights_w = (b&1) + 2;
            weights_h = a + 2;

        } else {
            // A2_B6
            weights_w = a + 2;
            weights_h = (b&1) + 6;
        }
    } else {
        // second 5 rows in table 11
        R |= _bits(2, 2) << 1;
        uint mode_bits = _bits(5, 4);
        if ( (mode_bits & 0xc) == 0 ) {
            // 12_A2
            die_assert( _bits(0, 4) != 0, "bad ASTC block mode" );
            weights_w = 12;
            weights_h = a + 2;
        } else if ( (mode_bits & 0xc) == 4 ) {
            // A2_12
            weights_w = a + 2;
            weights_h = 12;
        } else if ( mode_bits == 0xc ) {
            // 6_10
            weights_w = 6;
            weights_h = 10;
        } else if ( mode_bits == 0xd ) {
            // 10_6
            weights_w = 10;
            weights_h = 6;
        } else if ( (mode_bits & 0xc) == 8 ) {
            // A6_B6
            b = _bits(9, 2); 
            plane_cnt = 1;
            H = 0;
            weights_w = a + 6;
            weights_h = b + 6;
        } else {
            // illegal
            die_assert( false, "illegal ASTC block mode" );
            weights_w = 0;
            weights_h = 0;
        }
    }

    weights_qmode = (R - 2) + 6 * H;

    die_assert( plane_cnt == 1 || partition_cnt <= 3, "astc: dual-plane mode must have no more than 3 partitions" );
}

inline uint Model::Texture::astc_decode_partition( const unsigned char * bdata, uint s, uint t, uint r, uint partition_cnt, uint texel_cnt )
{
    //---------------------------------------------------------------
    // See section 3.15.
    //---------------------------------------------------------------
    if ( partition_cnt <= 1 ) return 0;

    uint seed = _bits(13, 10);

    if ( texel_cnt < 31 ) {
        s <<= 1;
        t <<= 1;
        r <<= 1;
    }

    seed += (partition_cnt - 1) * 1024;

    // hash52 inlined here:
    uint32_t rnum = seed;
    rnum ^= rnum >> 15;
    rnum -= rnum << 17;
    rnum += rnum << 7;
    rnum += rnum << 4;
    rnum ^= rnum >> 5;
    rnum += rnum << 16;
    rnum ^= rnum >> 7;
    rnum ^= rnum >> 3;
    rnum ^= rnum << 6;
    rnum ^= rnum >> 17;

    uint8_t seed1 = rnum & 0xF;
    uint8_t seed2 = (rnum >> 4) & 0xF;
    uint8_t seed3 = (rnum >> 8) & 0xF;
    uint8_t seed4 = (rnum >> 12) & 0xF;
    uint8_t seed5 = (rnum >> 16) & 0xF;
    uint8_t seed6 = (rnum >> 20) & 0xF;
    uint8_t seed7 = (rnum >> 24) & 0xF;
    uint8_t seed8 = (rnum >> 28) & 0xF;
    uint8_t seed9 = (rnum >> 18) & 0xF;
    uint8_t seed10 = (rnum >> 22) & 0xF;
    uint8_t seed11 = (rnum >> 26) & 0xF;
    uint8_t seed12 = ((rnum >> 30) | (rnum << 2)) & 0xF;

    seed1 *= seed1;
    seed2 *= seed2;
    seed3 *= seed3;
    seed4 *= seed4;
    seed5 *= seed5;
    seed6 *= seed6;
    seed7 *= seed7;
    seed8 *= seed8;
    seed9 *= seed9;
    seed10 *= seed10;
    seed11 *= seed11;
    seed12 *= seed12;

    uint sh1, sh2, sh3;
    if ( seed & 1 ) {
      sh1 = (seed & 2 ? 4 : 5);
      sh2 = (partition_cnt == 3 ? 6 : 5);
    } else {
      sh1 = (partition_cnt == 3 ? 6 : 5);
      sh2 = (seed & 2 ? 4 : 5);
    }
    sh3 = (seed & 0x10) ? sh1 : sh2;

    seed1 >>= sh1;
    seed2 >>= sh2;
    seed3 >>= sh1;
    seed4 >>= sh2;
    seed5 >>= sh1;
    seed6 >>= sh2;
    seed7 >>= sh1;
    seed8 >>= sh2;
    seed9 >>= sh3;
    seed10 >>= sh3;
    seed11 >>= sh3;
    seed12 >>= sh3;

    uint a = seed1 * s + seed2 * t + seed11 * r + (rnum >> 14);
    uint b = seed3 * s + seed4 * t + seed12 * r + (rnum >> 10);
    uint c = seed5 * s + seed6 * t + seed9  * r + (rnum >> 06);
    uint d = seed7 * s + seed8 * t + seed10 * r + (rnum >> 02);

    a &= 0x3F;
    b &= 0x3F;
    c &= 0x3F;
    d &= 0x3F;

    if ( partition_cnt <= 3 ) d = 0;
    if ( partition_cnt <= 2 ) c = 0;

    uint p;
    if ( a >= b && a >= c && a >= d ) {
        p = 0;
    } else if ( b >= c && b >= d ) {
        p = 1;
    } else if ( c >= d ) {
        p = 2;
    } else {
        p = 3;
    }

    mdout << "astc_decode_partition: partition=" << p << 
                        " rnum=0x" << std::hex << rnum << std::dec << " seed=" << seed << " sh1=" << sh1 << " sh2=" << sh2 << 
                        " seed1=" << int(seed1) << " seed2=" << int(seed2) << " seed3=" << int(seed3) << " seed4=" << int(seed4) <<
                        " seed5=" << int(seed5) << " seed6=" << int(seed6) << " seed7=" << int(seed7) << " seed8=" << int(seed8) <<
                        " pa=" << a << " pb=" << b << " pc=" << c << " pd=" << d << "\n";
    return p;
}

inline void Model::Texture::astc_decode_color_endpoint_mode( const unsigned char * bdata, uint plane_cnt, uint partition, uint partition_cnt, 
                                                             uint weights_bit_cnt, uint below_weights_start, 
                                                             uint& cem, uint& qmode, uint& trits, uint& quints, uint& bits, 
                                                             uint& endpoints_cnt, uint& endpoints_start, uint& endpoints_v_first, uint& endpoints_v_cnt )
{
    //---------------------------------------------------------------
    // See Section 3.5.
    // See Tables 13, 14, and 15.
    // This stuff is the the most confusing part of ASTC.
    // We follow the ARM astcenc physical_to_symbolic() code here.
    //
    // Determine the Color Endpoint Mode (CEM) for our partition.
    //---------------------------------------------------------------
    uint encoded_type_highpart_bit_cnt = 0;
    uint cems[4] = { 0, 0, 0, 0 };
    uint encoded_type = 0;
    uint base_class = 0;
    if ( partition_cnt == 1 ) {
        //---------------------------------------------------------------
        // One partition, so one set of CEM
        //---------------------------------------------------------------
        cems[0] = _bits(13, 4);
    } else {
        encoded_type_highpart_bit_cnt = (3 * partition_cnt) - 4;
        below_weights_start -= encoded_type_highpart_bit_cnt;
        //---------------------------------------------------------------
        // This is easier when we concatenate all the bits.
        //---------------------------------------------------------------
        encoded_type = _bits(23, 6) | (_bits(below_weights_start, encoded_type_highpart_bit_cnt) << 6);
        base_class = encoded_type & 3;
        if ( base_class == 0 ) {
            //---------------------------------------------------------------
            // All partitions have same CEM.
            //---------------------------------------------------------------
            for( uint i = 0; i < partition_cnt; i++ )
            {
                cems[i] = (encoded_type >> 2) & 0xf;
            }
            below_weights_start += encoded_type_highpart_bit_cnt;
        } else {
            //---------------------------------------------------------------
            // See Table 15.
            //---------------------------------------------------------------
            uint bit_pos = 2;
            for( uint i = 0; i < partition_cnt; i++, bit_pos++ )
            {
                cems[i] = (((encoded_type >> bit_pos) & 1) + base_class - 1) << 2;
            }
            for( uint i = 0; i < partition_cnt; i++, bit_pos += 2 )
            {
                cems[i] |= (encoded_type >> bit_pos) & 0x3;
            }
        }
    }
    cem = cems[partition];  

    //---------------------------------------------------------------
    // Determine number of integers we need to unpack if we were getting
    // all of the endpoint pairs.  
    //
    // Also, record the v_first and v_cnt for our partition.
    //
    // Then figure out the number of bits available for color endpoints.
    //
    // From that, we can then get the color quantization mode which 
    // is used by table lookups to unquantize color values.
    //---------------------------------------------------------------
    endpoints_cnt = 0;
    for( uint i = 0; i < partition_cnt; i++ )
    {
        uint v_cnt = 2 * (cems[i] >> 2) + 2;
        if ( i == partition ) {
            endpoints_v_first = endpoints_cnt;
            endpoints_v_cnt   = v_cnt;
        }
        endpoints_cnt += v_cnt;
    }
    die_assert( endpoints_cnt <= 18, "astc: too many color endpoint pairs" );

    int color_bits_cnt = (partition_cnt <= 1) ? (115 - 4) : (113 - 4 - 10);
    color_bits_cnt -= weights_bit_cnt + ((base_class != 0) ? encoded_type_highpart_bit_cnt : 0) + (plane_cnt-1)*2;
    if ( color_bits_cnt < 0 ) color_bits_cnt = 0;

    qmode = astc_color_quantization_mode[endpoints_cnt/2][color_bits_cnt];

    endpoints_start = (partition_cnt == 1) ? 17 : 29;
    mdout << "astc: cem=" << cem << 
                    " partition_cnt=" << partition_cnt << 
                    " weights_bit_cnt=" << weights_bit_cnt << 
                    " encoded_type_highpart_bit_cnt=" << encoded_type_highpart_bit_cnt << 
                    " plane_cnt=" << plane_cnt << 
                    " color_bits_cnt=" << color_bits_cnt << 
                    " cem_qmode=" << qmode << 
                    " endpoints_cnt=" << endpoints_cnt << 
                    " endpoints_start=" << endpoints_start << 
                    " endpoints_v_cnt=" << endpoints_v_cnt << 
                    " endpoints_v_first=" << endpoints_v_first << 
                    " encoded_type=" << encoded_type << 
                    " base_class=" << base_class << 
                    " cem0=" << cems[0] << " cem1=" << cems[1] << " cem2=" << cems[2] << " cem3=" << cems[3] << "\n";

    die_assert( qmode < 22, "astc: invalid color quantization mode, not < 22\n" );

    //---------------------------------------------------------------
    // Figure out the number of trits, quints, and bits for our partition's CEM.
    // Use table lookup.
    //---------------------------------------------------------------
    const ASTC_Range_Encoding& enc = astc_range_encodings[qmode];
    trits  = enc.trits;
    quints = enc.quints;
    bits   = enc.bits; 
    mdout << "astc: chosen range encoding: qmode=" << qmode << " max_p1=" << enc.max_p1 << " trits=" << trits << " quints=" << quints << " bits=" << bits << "\n";
}

inline void Model::Texture::astc_decode_color_endpoints( const unsigned char * bdata, uint cem, uint cem_qmode, 
                                                         uint trits, uint quints, uint bits, 
                                                         uint endpoints_cnt, uint endpoints_start, uint endpoints_v_first, uint endpoints_v_cnt,
                                                         uint endpoints[4][2] )
{
    //---------------------------------------------------------------
    // See Section 3.8 and table 20.
    // 
    // Read the unquantized v0, v1, etc. values starting at endpoints_v_first.
    // Keep them around as signed values.
    //---------------------------------------------------------------
    die_assert( endpoints_v_cnt <= 8, "astc: endpoints_v_cnt is > 8" );
    int  v[8];
    for( uint i = 0; i < endpoints_v_cnt; i++ )
    {
        //---------------------------------------------------------------
        // Read quantized value.
        // Unquantize it.
        //---------------------------------------------------------------
        uint ii = endpoints_v_first + i;
        uint qv = astc_decode_integer( bdata, endpoints_start, endpoints_cnt, ii, trits, quints, bits );
        v[i] = astc_color_unquantized[cem_qmode][qv];
        mdout << "astc_decode_color_endpoints: v" << i << "=" << v[i] << " v_first=" << endpoints_v_first << " v_cnt=" << endpoints_v_cnt << " i=" << i << " ii=" << ii << " qv=" << qv << "\n";
    }

    //---------------------------------------------------------------
    // Precompute these for debug.
    //---------------------------------------------------------------
    int v_bts[8] = { v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7] };
    astc_decode_bit_transfer_signed( v_bts[1], v_bts[0] );
    astc_decode_bit_transfer_signed( v_bts[3], v_bts[2] );
    astc_decode_bit_transfer_signed( v_bts[5], v_bts[4] );
    astc_decode_bit_transfer_signed( v_bts[7], v_bts[6] );

    for( uint i = 0; i < endpoints_v_cnt; i++ )
    {
        int bts = (v_bts[i] < 0) ? (0x200 + v_bts[i]) : v_bts[i];
        mdout << "astc_decode_color_endpoints: v" << i << "=" << v[i] << " v" << i << "_bts=" << bts << " real_bts=" << v_bts[i] << "\n";
    }

    //---------------------------------------------------------------
    // Now use those values according to the CEM.
    // Not supporting HDR right now.
    //---------------------------------------------------------------
    switch( cem )
    {
        case 0:
        {
            // Luminance
            endpoints[0][0] = v[0];
            endpoints[0][1] = v[1];
            endpoints[1][0] = v[0];
            endpoints[1][1] = v[1];
            endpoints[2][0] = v[0];
            endpoints[2][1] = v[1];
            endpoints[3][0] = 0xff;
            endpoints[3][1] = 0xff;
            break;
        }

        case 1:
        {
            // Luminance Delta
            uint L0 = (v[0] >> 2) | (v[1] & 0xc0);
            uint L1 = L0 + (v[1] & 0x3f);
            if ( L1 > 0xff ) L1 = 0xff;
            endpoints[0][0] = L0;
            endpoints[0][1] = L1;
            endpoints[1][0] = L0;
            endpoints[1][1] = L1;
            endpoints[2][0] = L0;
            endpoints[2][1] = L1;
            endpoints[3][0] = 0xff;
            endpoints[3][1] = 0xff;
            break;
        }
           
        case 4:
        {
            // Luminance Alpha
            endpoints[0][0] = v[0];
            endpoints[0][1] = v[1];
            endpoints[1][0] = v[0];
            endpoints[1][1] = v[1];
            endpoints[2][0] = v[0];
            endpoints[2][1] = v[1];
            endpoints[3][0] = v[2];
            endpoints[3][1] = v[3];
            break;
        }

        case 5:
        {
            // Luminance Alpha Delta
            endpoints[0][0] = astc_decode_clamp_unorm8( v_bts[0] );
            endpoints[0][1] = astc_decode_clamp_unorm8( v_bts[0] + v_bts[1] );
            endpoints[1][0] = endpoints[0][0];
            endpoints[1][1] = endpoints[0][1];
            endpoints[2][0] = endpoints[0][0];
            endpoints[2][1] = endpoints[0][1];
            endpoints[3][0] = astc_decode_clamp_unorm8( v_bts[2] );
            endpoints[3][1] = astc_decode_clamp_unorm8( v_bts[2] + v_bts[3] );
            break;
        }

        case 6:
        {
            // RGB Scale
            endpoints[0][0] = (v[0]*v[3]) >> 8;
            endpoints[0][1] = v[0];
            endpoints[1][0] = (v[1]*v[3]) >> 8;
            endpoints[1][1] = v[1];
            endpoints[2][0] = (v[2]*v[3]) >> 8;
            endpoints[2][1] = v[2];
            endpoints[3][0] = 0xff;
            endpoints[3][1] = 0xff;
            break;
        }

        case 8:
        {
            // RGB 
            int s0 = v[0] + v[2] + v[4];
            int s1 = v[1] + v[3] + v[5];
            mdout << "astc_decode_color_endpoints: RGB: rgb0=[" << v[0] << "," << v[2] << "," << v[4] << 
                                                     "] rgb1=[" << v[1] << "," << v[3] << "," << v[5] << "]\n";
            if ( s1 >= s0 ) {
                // no blue contraction
                endpoints[0][0] = astc_decode_clamp_unorm8( v[0] );
                endpoints[0][1] = astc_decode_clamp_unorm8( v[1] );
                endpoints[1][0] = astc_decode_clamp_unorm8( v[2] );
                endpoints[1][1] = astc_decode_clamp_unorm8( v[3] );
                endpoints[2][0] = astc_decode_clamp_unorm8( v[4] );
                endpoints[2][1] = astc_decode_clamp_unorm8( v[5] );
                endpoints[3][0] = 0xff;
                endpoints[3][1] = 0xff;
                mdout << "astc_decode_color_endpoints: RGB: no blue contraction\n";
            } else {
                // blue contraction
                int a = 0xff;
                astc_decode_blue_contract( v[1], v[3], v[5], a );
                endpoints[0][0] = astc_decode_clamp_unorm8( v[1] );
                endpoints[1][0] = astc_decode_clamp_unorm8( v[3] );
                endpoints[2][0] = astc_decode_clamp_unorm8( v[5] );
                endpoints[3][0] = a;

                a = 0xff;
                astc_decode_blue_contract( v[0], v[2], v[4], a );
                endpoints[0][1] = astc_decode_clamp_unorm8( v[0] );
                endpoints[1][1] = astc_decode_clamp_unorm8( v[2] );
                endpoints[2][1] = astc_decode_clamp_unorm8( v[4] );
                endpoints[3][1] = a;

                mdout << "astc_decode_color_endpoints: RGB: blue contraction: rgb0=[" << v[0] << "," << v[2] << "," << v[4] << 
                                                                           "] rgb1=[" << v[1] << "," << v[3] << "," << v[5] << "]\n";
            }
            break;
        }

        case 9:
        {
            // RGB Delta
            int r = v_bts[0] + v_bts[1];
            int g = v_bts[2] + v_bts[3];
            int b = v_bts[4] + v_bts[5];
            int a = 0xff;
            endpoints[3][0] = a;
            endpoints[3][1] = a;
            int rgbsum = v_bts[1] + v_bts[3] + v_bts[5];
            mdout << "astc_decode_color_endpoints: RGB Delta: brgb0=[" << v_bts[0] << "," << v_bts[2] << "," << v_bts[4] << "]" <<
                                                            " brgb1=[" << r    << "," << g    << "," << b    << "] rgbsum=" << rgbsum << "\n";
            if ( rgbsum >= 0 ) {
                endpoints[0][0] = astc_decode_clamp_unorm8( v_bts[0] );
                endpoints[1][0] = astc_decode_clamp_unorm8( v_bts[2] );
                endpoints[2][0] = astc_decode_clamp_unorm8( v_bts[4] );
                endpoints[0][1] = astc_decode_clamp_unorm8( r );
                endpoints[1][1] = astc_decode_clamp_unorm8( g );
                endpoints[2][1] = astc_decode_clamp_unorm8( b );
            } else {
                astc_decode_blue_contract( r, g, b, a );
                endpoints[0][0] = astc_decode_clamp_unorm8( r );
                endpoints[1][0] = astc_decode_clamp_unorm8( g );
                endpoints[2][0] = astc_decode_clamp_unorm8( b );

                r = v_bts[0];
                g = v_bts[2];
                b = v_bts[4];
                astc_decode_blue_contract( r, g, b, a );
                endpoints[0][1] = astc_decode_clamp_unorm8( r );
                endpoints[1][1] = astc_decode_clamp_unorm8( g );
                endpoints[2][1] = astc_decode_clamp_unorm8( b );
            }
            break;
        }

        case 10:
        {
            // RGB Scale Alpha
            endpoints[0][0] = (v[0]*v[3]) >> 8;
            endpoints[0][1] = v[0];
            endpoints[1][0] = (v[1]*v[3]) >> 8;
            endpoints[1][1] = v[1];
            endpoints[2][0] = (v[2]*v[3]) >> 8;
            endpoints[2][1] = v[2];
            endpoints[3][0] = v[4];
            endpoints[3][1] = v[5];
            break;
        }

        case 12:
        {
            // RGBA
            if ( (v[1] + v[3] + v[5]) >= (v[0] + v[2] + v[4]) ) {
                endpoints[0][0] = v[0];
                endpoints[1][0] = v[2];
                endpoints[2][0] = v[4];
                endpoints[3][0] = v[6];

                endpoints[0][1] = v[1];
                endpoints[1][1] = v[3];
                endpoints[2][1] = v[5];
                endpoints[3][1] = v[7];
            } else {
                astc_decode_blue_contract( v[1], v[3], v[5], v[7] );
                endpoints[0][0] = v[1];
                endpoints[1][0] = v[3];
                endpoints[2][0] = v[5];
                endpoints[3][0] = v[7];

                astc_decode_blue_contract( v[0], v[2], v[4], v[6] );
                endpoints[0][1] = v[0];
                endpoints[1][1] = v[2];
                endpoints[2][1] = v[4];
                endpoints[3][1] = v[6];
            }
            break;
        }

        case 13:
        {
            // RGBA Delta
            if ( (v_bts[1] + v_bts[3] + v_bts[5]) >= 0 ) {
                endpoints[0][0] = v_bts[0];
                endpoints[1][0] = v_bts[2];
                endpoints[2][0] = v_bts[4];
                endpoints[3][0] = v_bts[6];

                endpoints[0][1] = astc_decode_clamp_unorm8( v_bts[0] + v_bts[1] );
                endpoints[1][1] = astc_decode_clamp_unorm8( v_bts[2] + v_bts[3] );
                endpoints[2][1] = astc_decode_clamp_unorm8( v_bts[4] + v_bts[5] );
                endpoints[3][1] = astc_decode_clamp_unorm8( v_bts[6] + v_bts[7] );
            } else {
                int r = v_bts[0] + v_bts[1];
                int g = v_bts[2] + v_bts[3];
                int b = v_bts[4] + v_bts[5];
                int a = v_bts[6] + v_bts[7];
                astc_decode_blue_contract( r, g, b, a );
                endpoints[0][0] = r;
                endpoints[1][0] = g;
                endpoints[2][0] = b;
                endpoints[3][0] = a;

                astc_decode_blue_contract( v_bts[0], v_bts[2], v_bts[4], v_bts[6] );
                endpoints[0][1] = v_bts[0];
                endpoints[1][1] = v_bts[2];
                endpoints[2][1] = v_bts[4];
                endpoints[3][1] = v_bts[6];
            }
            break;
        }

        default:
        {
#ifdef ASTC_HACK_HDR
            // temporary hack so programs don't die, will go away soon
            endpoints[0][0] = v[0];
            endpoints[1][0] = v[2];
            endpoints[2][0] = v[4];
            endpoints[3][0] = v[6];

            endpoints[0][1] = v[1];
            endpoints[1][1] = v[3];
            endpoints[2][1] = v[5];
            endpoints[3][1] = v[7];
#else
            die_assert( false, "ASTC HDR color endpoint modes are not yet supported" );
#endif
            break;
        }
    }

    mdout << "astc_decode_color_endpoints: final: rgba0=[" << endpoints[0][0] << "," << endpoints[1][0] << "," << endpoints[2][0] << "," << endpoints[3][0] << 
                                               "] rgba1=[" << endpoints[0][1] << "," << endpoints[1][1] << "," << endpoints[2][1] << "," << endpoints[3][1] << "]\n";
}

inline void Model::Texture::astc_decode_bit_transfer_signed( int& a, int& b )
{
    //---------------------------------------------------------------
    // See code below Table 20.
    //---------------------------------------------------------------
    uint a_orig = a;
    uint b_orig = b;
    b |= (a & 0x80) << 1;
    a &= 0x7f;
    if ( (a & 0x40) != 0 ) a -= 0x80;
    b >>= 1;
    a >>= 1;
    mdout << "bts: a_orig=" << a_orig << " b_orig=" << b_orig << " a=" << a << " b=" << b << "\n";
}

inline void Model::Texture::astc_decode_blue_contract( int& r, int& g, int& b, int& a )
{
    //---------------------------------------------------------------
    // See code below Table 20.
    // b and a don't change.
    //---------------------------------------------------------------
    r = (r+b) >> 1;
    g = (g+b) >> 1;
    (void)a;
}

inline uint Model::Texture::astc_decode_clamp_unorm8( int v )
{
    //---------------------------------------------------------------
    // Clamp to 0 .. 255 range.
    //---------------------------------------------------------------
    return (v < 0) ? 0 : (v > 255) ? 255 : v;
}

inline void Model::Texture::astc_decode_weights( const unsigned char * rbdata, uint weights_start, uint plane_cnt,
                                                 uint blk_w, uint blk_h, uint blk_d, uint weights_w, uint weights_h, uint weights_d,
                                                 uint s, uint t, uint r, uint trits, uint quints, uint bits, uint qmode,
                                                 uint plane_weights[2] )
{

    die_assert( blk_d == 1 && weights_d == 1 && r == 0, "can't handle 3D textures right now" );
    mdout << "astc_decode_weights: blk_whd=[" << blk_w << "," << blk_h << "," << blk_d << "]" <<
                                 " weights_whd=[" << weights_w << "," << weights_h << "," << weights_d << "]\n";

    //---------------------------------------------------------------
    // DECIMATE
    //
    // Look up the decimation encodings for this combination of 
    // block and weights dimensions.  It will tell us which
    // weights to pull out for this texel and decimation factors.
    //---------------------------------------------------------------
    uint texel_i = s + t*blk_w;
    const ASTC_Block_Weights_Decimation_Encoding& decimations = astc_block_weights_decimation_encoding( blk_w, blk_h, weights_w, weights_h );
    const ASTC_Texel_Decimation_Encoding&         decimation  = decimations[texel_i];
    uint weight_cnt = decimation.weight_cnt;
    uint real_weight_cnt = plane_cnt * weights_w * weights_h * weights_d;
    mdout << "astc_decode_weights: texel_i=" << texel_i << " decimation_weight_cnt=" << weight_cnt << " real_weight_cnt=" << real_weight_cnt << 
                                 " weights_start=" << weights_start << "\n";
    for( uint p = 0; p < plane_cnt; p++ )
    {
        uint sum = 8;
        for( uint i = 0; i < weight_cnt; i++ )
        {
            //---------------------------------------------------------------
            // Read quantized weight wi.
            // Then unquantize it using a table lookup.
            //---------------------------------------------------------------
            uint wi  = decimation.weight_i[i];
            uint wf  = decimation.weight_factor[i];
            uint wii = plane_cnt*wi + p;
            uint qw  = astc_decode_integer( rbdata, weights_start, real_weight_cnt, wii, trits, quints, bits );
            die_assert( qw < 32, "astc_decode_weights: qw=" + std::to_string(qw) + " which is >= 32" );
            uint w   = astc_weight_unquantized[qmode][qw];
            sum += w * wf;
            mdout << "astc_decode_weights: p=" << p << " i=" << i << " wi=" << wi << " wf=" << wf << " wii=" << wii << 
                                        " qmode=" << qmode << " qw=" << qw << " w=" << w << " sum=" << sum << "\n";
        }
        plane_weights[p] = sum >> 4;
        mdout << "astc_decode_weights: plane_weights[" << p << "]=" << plane_weights[p] << "\n";
    }
}

const Model::Texture::ASTC_Block_Weights_Decimation_Encoding& 
                    Model::Texture::astc_block_weights_decimation_encoding( uint blk_w, uint blk_h, uint weights_w, uint weights_h )
{
    //---------------------------------------------------------------
    // If the table for this combo has already been created, then 
    // return it quickly.
    //---------------------------------------------------------------
    uint blk_i     = 16*blk_h     + blk_w;
    uint weights_i = 16*weights_h + weights_w;
    uint texel_cnt = blk_w * blk_h;
    if ( astc_weight_decimation_encodings.size() > blk_i ) {
        ASTC_Block_Decimation_Encoding& blk_encodings = astc_weight_decimation_encodings[blk_i];
        if ( blk_encodings.size() > weights_i ) {
            Model::Texture::ASTC_Block_Weights_Decimation_Encoding& encodings = blk_encodings[weights_i];
            if ( encodings.size() != 0 ) {
                die_assert( encodings.size() == texel_cnt, "encodings inconsistency detected" );
                return encodings;
            }
        }
    }

    //---------------------------------------------------------------
    // No luck.  
    // First make sure all the parent arrays are there.
    //---------------------------------------------------------------
    if ( astc_weight_decimation_encodings.size() <= blk_i ) {
        astc_weight_decimation_encodings.resize( blk_i+1 );
    }
    ASTC_Block_Decimation_Encoding& blk_encodings = astc_weight_decimation_encodings[blk_i];
    if ( blk_encodings.size() <= weights_i ) {
        blk_encodings.resize( weights_i+1 );
    }
    Model::Texture::ASTC_Block_Weights_Decimation_Encoding& encodings = blk_encodings[weights_i];
    
    //---------------------------------------------------------------
    // Now fill in this combo.
    //---------------------------------------------------------------
    mdout << "Decimations for blk_whd=[" << blk_w << "," << blk_h << ",1]" <<
             " weights_whd=[" << weights_w << "," << weights_h << ",1] weights_i=" << weights_i << "\n";
    die_assert( encodings.size() == 0, "encodings should have been empty" );
    encodings.resize( texel_cnt );

    for( uint y = 0; y < blk_h; y++ )
    {
        for( uint x = 0; x < blk_w; x++ )
        {
            uint texel_i = y*blk_w + x;

            uint x_weight = (((1024 + blk_w / 2) / (blk_w - 1)) * x * (weights_w - 1) + 32) >> 6;
            uint y_weight = (((1024 + blk_h / 2) / (blk_h - 1)) * y * (weights_h - 1) + 32) >> 6;

            uint x_weight_frac = x_weight & 0xF;
            uint y_weight_frac = y_weight & 0xF;
            uint x_weight_int  = x_weight >> 4;
            uint y_weight_int  = y_weight >> 4;

            uint qweight[4];
            qweight[0] = x_weight_int + y_weight_int * weights_w;
            qweight[1] = qweight[0] + 1;
            qweight[2] = qweight[0] + x_weight;
            qweight[3] = qweight[2] + 1;

            uint prod = x_weight_frac * y_weight_frac;

            uint weight[4];
            weight[3] = (prod + 8) >> 4;
            weight[1] = x_weight_frac - weight[3];
            weight[2] = y_weight_frac - weight[3];
            weight[0] = 16 - x_weight_frac - y_weight_frac + weight[3];

            mdout << "    [" << x << "," << y << "]: x_weight=" << x_weight << " y_weight=" << y_weight << " prod=" << prod << 
                        " qweight[]=[" << qweight[0] << "," << qweight[1] << "," << qweight[2] << "," << qweight[3] << "]" << 
                        " weight[]=["  << weight[0]  << "," << weight[1]  << "," << weight[2]  << "," << weight[3]  << "]\n";

            ASTC_Texel_Decimation_Encoding& encoding = encodings[texel_i];
            encoding.weight_cnt = 0;
            for( uint i = 0; i < 4; i++ )
            {
                if ( weight[i] != 0 ) {
                    encoding.weight_i[encoding.weight_cnt] = qweight[i];
                    encoding.weight_factor[encoding.weight_cnt++] = weight[i];
                }
            }
        }
    }
    return encodings;
}

inline uint Model::Texture::astc_decode_integer( const unsigned char * bdata, uint start, uint cnt, uint vi, uint trits, uint quints, uint bits )
{
    //---------------------------------------------------------------
    // See SPEC Section 3.6.
    // See decode_ise() from open-source ASTC implementation.
    // This is tricky stuff so we follow decode_ise() rather than the spec.
    //
    // Note: bits == n == number of LSBs in each value
    //       trits and quints can add other MSBs.
    //
    // One difference between this routine and decode_ise() is that
    // we care about only one integer, so we skip all values until value vi.
    // We currently use loops to do this skipping, but later we will likely 
    // add lookup tables (LUTs) to avoid these loops.
    //---------------------------------------------------------------
    mdout << "astc_decode_integer: start=" << start << " cnt=" << cnt << " vi=" << vi << 
             " trits=" << trits << " quints=" << quints << " bits=" << bits << " dat=" << astc_dat_str( bdata ) << "\n";
    const uint TQ_SIZE = 22;
    uint tq_blocks[TQ_SIZE];		// trit-blocks or quint-blocks
    for( uint i = 0; i < TQ_SIZE; i++ )
        tq_blocks[i] = 0;

    uint lcounter = 0;
    uint hcounter = 0;

    // collect bits for element vi, as well as bits for any trit-blocks and quint-blocks.
    uint v;
#ifdef ASTC_HACK_HDR
    vi = vi % cnt;
#else
    die_assert( vi < cnt, "astc_decode_integer vi >= cnt" );
#endif
    for( uint i = 0; i < cnt; i++ )  // should be able to short-circuit this better later
    {
        if ( i == vi ) {
            v = _astc_bits( bdata, start, bits );
            mdout << "astc_decode_integer: i=" << i << " initial v=" << v << " start=" << start << " bits=" << bits << "\n";
        }
        start += bits;

        if( trits ) {
            static const uint bits_to_read[5]  = { 2, 2, 1, 2, 1 };
            static const uint block_shift[5]   = { 0, 2, 4, 5, 7 };
            static const uint next_lcounter[5] = { 1, 2, 3, 4, 0 };
            static const uint hcounter_incr[5] = { 0, 0, 0, 0, 1 };
            uint t = _astc_bits( bdata, start, bits_to_read[lcounter] );
            die_assert( hcounter < TQ_SIZE, "trits tq_blocks[] overflow vi=" + std::to_string(vi) );
            tq_blocks[hcounter] |= t << block_shift[lcounter];
            mdout << "astc_decode_integer: i=" << i << " trit t=" << t << " start=" << start << " bits_to_read=" << bits_to_read[lcounter] <<
                                         " hcounter=" << hcounter << " lcounter=" << lcounter << 
                                         " tq_blocks[" << hcounter << "]=" << tq_blocks[hcounter] << "\n";
            start += bits_to_read[lcounter];
            hcounter += hcounter_incr[lcounter];
            lcounter  = next_lcounter[lcounter];
        }

        if( quints ) {
            static const uint bits_to_read[3]  = { 3, 2, 2 };
            static const uint block_shift[3]   = { 0, 3, 5 };
            static const uint next_lcounter[3] = { 1, 2, 0 };
            static const uint hcounter_incr[3] = { 0, 0, 1 };
            uint q = _astc_bits( bdata, start, bits_to_read[lcounter] );
            die_assert( hcounter < TQ_SIZE, "quints tq_blocks[] overflow vi=" + std::to_string(vi) );
            tq_blocks[hcounter] |= q << block_shift[lcounter];
            mdout << "astc_decode_integer: i=" << i << " quint q=" << q << " start=" << start << " bits_to_read=" << bits_to_read[lcounter] <<
                                         " hcounter=" << hcounter << " lcounter=" << lcounter << 
                                         " tq_block[" << hcounter << "]=" << tq_blocks[hcounter] << "\n";
            start += bits_to_read[lcounter];
            hcounter += hcounter_incr[lcounter];
            lcounter  = next_lcounter[lcounter];
        }
    }

    // unpack trit or quint block that we need for value vi
    if( trits )
    {
        uint ti = vi / 5;
        uint to = vi % 5;
        die_assert( ti < TQ_SIZE, "reading beyond end of trits tq_blocks[] ti=" + std::to_string(ti) );
        uint tq = tq_blocks[ti];
        uint tt = astc_integer_trits[tq][to];
        v |= tt << bits;
        mdout << "astc_decode_integer: trits[" << tq << "][" << to << "]=tt=" << tt << " << " << bits << " = " << v << "\n";
    } else if ( quints ) {
        uint qi = vi / 3;
        uint qo = vi % 3;
        die_assert( qi < TQ_SIZE, "reading beyond end of quints tq_blocks[] qi=" + std::to_string(qi) );
        uint tq = tq_blocks[qi];
        uint qq = astc_integer_quints[tq][qo];
        v |= qq << bits;
        mdout << "astc_decode_integer: quints[" << tq << "][" << qo << "]=qq=" << qq << " << " << bits << " = " << v << "\n";
    } 

    return v;
}

bool Model::Texture::ASTC_Header::magic_is_good( void ) const
{
    uint64_t full_magic = (magic[0] << 0) | (magic[1] << 8) | (magic[2] << 16) | (magic[3] << 24);
    return full_magic == 0x5CA1AB13;
}

uint Model::Texture::ASTC_Header::size_x( void ) const
{
    return (xsize[0] << 0) | (xsize[1] << 8) | (xsize[2] << 16);
}

uint Model::Texture::ASTC_Header::size_y( void ) const
{
    return (ysize[0] << 0) | (ysize[1] << 8) | (ysize[2] << 16);
}

uint Model::Texture::ASTC_Header::size_z( void ) const
{
    return (zsize[0] << 0) | (zsize[1] << 8) | (zsize[2] << 16);
}

uint Model::Texture::ASTC_Header::blk_cnt_x( void ) const
{
    return (size_x() + blockdim_x - 1) / blockdim_x;
}

uint Model::Texture::ASTC_Header::blk_cnt_y( void ) const
{
    return (size_y() + blockdim_y - 1) / blockdim_y;
}

uint Model::Texture::ASTC_Header::blk_cnt_z( void ) const
{
    return (size_z() + blockdim_z - 1) / blockdim_z;
}

uint Model::Texture::ASTC_Header::blk_cnt( void ) const
{
    return 1 + blk_cnt_x()*blk_cnt_y()*blk_cnt_z();
}

uint Model::Texture::ASTC_Header::byte_cnt( void ) const
{
    return 16*blk_cnt();
}

/*
   This software is available under 2 licenses -- choose whichever you prefer.
   ------------------------------------------------------------------------------
   ALTERNATIVE A - MIT License
   Copyright (c) 2017 Sean Barrett
   Permission is hereby granted, free of charge, to any person obtaining a copy of
   this software and associated documentation files (the "Software"), to deal in
   the Software without restriction, including without limitation the rights to
   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is furnished to do
   so, subject to the following conditions:
   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   ------------------------------------------------------------------------------
   ALTERNATIVE B - Public Domain (www.unlicense.org)
   This is free and unencumbered software released into the public domain.
   Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
   software, either in source code form or as a compiled binary, for any purpose,
   commercial or non-commercial, and by any means.
   In jurisdictions that recognize copyright laws, the author or authors of this
   software dedicate any and all copyright interest in the software to the public
   domain. We make this dedication for the benefit of the public at large and to
   the detriment of our heirs and successors. We intend this dedication to be an
   overt act of relinquishment in perpetuity of all present and future rights to
   this software under copyright law.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

   stb_image.h - v2.06 - public domain image loader - http://nothings.org/stb_image.h
                       no warranty implied; use at your own risk

   You can #define STBI_ASSERT(x) before the #include to avoid using assert.h.
   And #define STBI_MALLOC, STBI_REALLOC, and STBI_FREE to avoid using malloc,realloc,free


   QUICK NOTES:
      Primarily of interest to game developers and other people who can
          avoid problematic images and only need the trivial interface

      JPEG baseline & progressive (12 bpc/arithmetic not supported, same as stock IJG lib)
      PNG 1/2/4/8-bit-per-channel (16 bpc not supported)

      TGA (not sure what subset, if a subset)
      BMP non-1bpp, non-RLE
      PSD (composited view only, no extra channels)

      GIF (*comp always reports as 4-channel)
      HDR (radiance rgbE format)
      PIC (Softimage PIC)
      PNM (PPM and PGM binary only)

      - decode from memory or through FILE (define STBI_NO_STDIO to remove code)
      - decode from arbitrary I/O callbacks
      - SIMD acceleration on x86/x64 (SSE2) and ARM (NEON)

   Full documentation under "DOCUMENTATION" below.


   Revision 2.00 release notes:

      - Progressive JPEG is now supported.

      - PPM and PGM binary formats are now supported, thanks to Ken Miller.

      - x86 platforms now make use of SSE2 SIMD instructions for
        JPEG decoding, and ARM platforms can use NEON SIMD if requested.
        This work was done by Fabian "ryg" Giesen. SSE2 is used by
        default, but NEON must be enabled explicitly; see docs.

        With other JPEG optimizations included in this version, we see
        2x speedup on a JPEG on an x86 machine, and a 1.5x speedup
        on a JPEG on an ARM machine, relative to previous versions of this
        library. The same results will not obtain for all JPGs and for all
        x86/ARM machines. (Note that progressive JPEGs are significantly
        slower to decode than regular JPEGs.) This doesn't mean that this
        is the fastest JPEG decoder in the land; rather, it brings it
        closer to parity with standard libraries. If you want the fastest
        decode, look elsewhere. (See "Philosophy" section of docs below.)

        See final bullet items below for more info on SIMD.

      - Added STBI_MALLOC, STBI_REALLOC, and STBI_FREE macros for replacing
        the memory allocator. Unlike other STBI libraries, these macros don't
        support a context parameter, so if you need to pass a context in to
        the allocator, you'll have to store it in a global or a thread-local
        variable.

      - Split existing STBI_NO_HDR flag into two flags, STBI_NO_HDR and
        STBI_NO_LINEAR.
            STBI_NO_HDR:     suppress implementation of .hdr reader format
            STBI_NO_LINEAR:  suppress high-dynamic-range light-linear float API

      - You can suppress implementation of any of the decoders to reduce
        your code footprint by #defining one or more of the following
        symbols before creating the implementation.

            STBI_NO_JPEG
            STBI_NO_PNG
            STBI_NO_BMP
            STBI_NO_PSD
            STBI_NO_TGA
            STBI_NO_GIF
            STBI_NO_HDR
            STBI_NO_PIC
            STBI_NO_PNM   (.ppm and .pgm)

      - You can request *only* certain decoders and suppress all other ones
        (this will be more forward-compatible, as addition of new decoders
        doesn't require you to disable them explicitly):

            STBI_ONLY_JPEG
            STBI_ONLY_PNG
            STBI_ONLY_BMP
            STBI_ONLY_PSD
            STBI_ONLY_TGA
            STBI_ONLY_GIF
            STBI_ONLY_HDR
            STBI_ONLY_PIC
            STBI_ONLY_PNM   (.ppm and .pgm)

         Note that you can define multiples of these, and you will get all
         of them ("only x" and "only y" is interpreted to mean "only x&y").

       - If you use STBI_NO_PNG (or _ONLY_ without PNG), and you still
         want the zlib decoder to be available, #define STBI_SUPPORT_ZLIB

      - Compilation of all SIMD code can be suppressed with
            #define STBI_NO_SIMD
        It should not be necessary to disable SIMD unless you have issues
        compiling (e.g. using an x86 compiler which doesn't support SSE
        intrinsics or that doesn't support the method used to detect
        SSE2 support at run-time), and even those can be reported as
        bugs so I can refine the built-in compile-time checking to be
        smarter.

      - The old STBI_SIMD system which allowed installing a user-defined
        IDCT etc. has been removed. If you need this, don't upgrade. My
        assumption is that almost nobody was doing this, and those who
        were will find the built-in SIMD more satisfactory anyway.

      - RGB values computed for JPEG images are slightly different from
        previous versions of stb_image. (This is due to using less
        integer precision in SIMD.) The C code has been adjusted so
        that the same RGB values will be computed regardless of whether
        SIMD support is available, so your app should always produce
        consistent results. But these results are slightly different from
        previous versions. (Specifically, about 3% of available YCbCr values
        will compute different RGB results from pre-1.49 versions by +-1;
        most of the deviating values are one smaller in the G channel.)

      - If you must produce consistent results with previous versions of
        stb_image, #define STBI_JPEG_OLD and you will get the same results
        you used to; however, you will not get the SIMD speedups for
        the YCbCr-to-RGB conversion step (although you should still see
        significant JPEG speedup from the other changes).

        Please note that STBI_JPEG_OLD is a temporary feature; it will be
        removed in future versions of the library. It is only intended for
        near-term back-compatibility use.


   Latest revision history:
      2.06  (2015-04-19) fix bug where PSD returns wrong '*comp' value
      2.05  (2015-04-19) fix bug in progressive JPEG handling, fix warning
      2.04  (2015-04-15) try to re-enable SIMD on MinGW 64-bit
      2.03  (2015-04-12) additional corruption checking
                         stbi_set_flip_vertically_on_load
                         fix NEON support; fix mingw support
      2.02  (2015-01-19) fix incorrect assert, fix warning
      2.01  (2015-01-17) fix various warnings
      2.00b (2014-12-25) fix STBI_MALLOC in progressive JPEG
      2.00  (2014-12-25) optimize JPEG, including x86 SSE2 & ARM NEON SIMD
                         progressive JPEG
                         PGM/PPM support
                         STBI_MALLOC,STBI_REALLOC,STBI_FREE
                         STBI_NO_*, STBI_ONLY_*
                         GIF bugfix
      1.48  (2014-12-14) fix incorrectly-named assert()
      1.47  (2014-12-14) 1/2/4-bit PNG support (both grayscale and paletted)
                         optimize PNG
                         fix bug in interlaced PNG with user-specified channel count

   See end of file for full revision history.


 ============================    Contributors    =========================

 Image formats                                Bug fixes & warning fixes
    Sean Barrett (jpeg, png, bmp)                Marc LeBlanc
    Nicolas Schulz (hdr, psd)                    Christpher Lloyd
    Jonathan Dummer (tga)                        Dave Moore
    Jean-Marc Lienher (gif)                      Won Chun
    Tom Seddon (pic)                             the Horde3D community
    Thatcher Ulrich (psd)                        Janez Zemva
    Ken Miller (pgm, ppm)                        Jonathan Blow
                                                 Laurent Gomila
                                                 Aruelien Pocheville
 Extensions, features                            Ryamond Barbiero
    Jetro Lauha (stbi_info)                      David Woo
    Martin "SpartanJ" Golini (stbi_info)         Martin Golini
    James "moose2000" Brown (iPhone PNG)         Roy Eltham
    Ben "Disch" Wenger (io callbacks)            Luke Graham
    Omar Cornut (1/2/4-bit PNG)                  Thomas Ruf
    Nicolas Guillemot (vertical flip)            John Bartholomew
                                                 Ken Hamada
 Optimizations & bugfixes                        Cort Stratton
    Fabian "ryg" Giesen                          Blazej Dariusz Roszkowski
    Arseny Kapoulkine                            Thibault Reuille
                                                 Paul Du Bois
                                                 Guillaume George
  If your name should be here but                Jerry Jansson
  isn't, let Sean know.                          Hayaki Saito
                                                 Johan Duparc
                                                 Ronny Chevalier
                                                 Michal Cichon
                                                 Tero Hanninen
                                                 Sergio Gonzalez
                                                 Cass Everitt
                                                 Engin Manap
                                                 Martins Mozeiko
                                                 Joseph Thomson
                                                 Phil Jordan

License:
   This software is in the public domain. Where that dedication is not
   recognized, you are granted a perpetual, irrevocable license to copy
   and modify this file however you want.

*/

// DOCUMENTATION
//
// Limitations:
//    - no 16-bit-per-channel PNG
//    - no 12-bit-per-channel JPEG
//    - no JPEGs with arithmetic coding
//    - no 1-bit BMP
//    - GIF always returns *comp=4
//
// Basic usage (see HDR discussion below for HDR usage):
//    int x,y,n;
//    unsigned char *data = stbi_load(filename, &x, &y, &n, 0);
//    // ... process data if not NULL ...
//    // ... x = width, y = height, n = # 8-bit components per pixel ...
//    // ... replace '0' with '1'..'4' to force that many components per pixel
//    // ... but 'n' will always be the number that it would have been if you said 0
//    stbi_image_free(data)
//
// Standard parameters:
//    int *x       -- outputs image width in pixels
//    int *y       -- outputs image height in pixels
//    int *comp    -- outputs # of image components in image file
//    int req_comp -- if non-zero, # of image components requested in result
//
// The return value from an image loader is an 'unsigned char *' which points
// to the pixel data, or NULL on an allocation failure or if the image is
// corrupt or invalid. The pixel data consists of *y scanlines of *x pixels,
// with each pixel consisting of N interleaved 8-bit components; the first
// pixel pointed to is top-left-most in the image. There is no padding between
// image scanlines or between pixels, regardless of format. The number of
// components N is 'req_comp' if req_comp is non-zero, or *comp otherwise.
// If req_comp is non-zero, *comp has the number of components that _would_
// have been output otherwise. E.g. if you set req_comp to 4, you will always
// get RGBA output, but you can check *comp to see if it's trivially opaque
// because e.g. there were only 3 channels in the source image.
//
// An output image with N components has the following components interleaved
// in this order in each pixel:
//
//     N=#comp     components
//       1           grey
//       2           grey, alpha
//       3           red, green, blue
//       4           red, green, blue, alpha
//
// If image loading fails for any reason, the return value will be NULL,
// and *x, *y, *comp will be unchanged. The function stbi_failure_reason()
// can be queried for an extremely brief, end-user unfriendly explanation
// of why the load failed. Define STBI_NO_FAILURE_STRINGS to avoid
// compiling these strings at all, and STBI_FAILURE_USERMSG to get slightly
// more user-friendly ones.
//
// Paletted PNG, BMP, GIF, and PIC images are automatically depalettized.
//
// ===========================================================================
//
// Philosophy
//
// stb libraries are designed with the following priorities:
//
//    1. easy to use
//    2. easy to maintain
//    3. good performance
//
// Sometimes I let "good performance" creep up in priority over "easy to maintain",
// and for best performance I may provide less-easy-to-use APIs that give higher
// performance, in addition to the easy to use ones. Nevertheless, it's important
// to keep in mind that from the standpoint of you, a client of this library,
// all you care about is #1 and #3, and stb libraries do not emphasize #3 above all.
//
// Some secondary priorities arise directly from the first two, some of which
// make more explicit reasons why performance can't be emphasized.
//
//    - Portable ("ease of use")
//    - Small footprint ("easy to maintain")
//    - No dependencies ("ease of use")
//
// ===========================================================================
//
// I/O callbacks
//
// I/O callbacks allow you to read from arbitrary sources, like packaged
// files or some other source. Data read from callbacks are processed
// through a small internal buffer (currently 128 bytes) to try to reduce
// overhead.
//
// The three functions you must define are "read" (reads some bytes of data),
// "skip" (skips some bytes of data), "eof" (reports if the stream is at the end).
//
// ===========================================================================
//
// SIMD support
//
// The JPEG decoder will try to automatically use SIMD kernels on x86 when
// supported by the compiler. For ARM Neon support, you must explicitly
// request it.
//
// (The old do-it-yourself SIMD API is no longer supported in the current
// code.)
//
// On x86, SSE2 will automatically be used when available based on a run-time
// test; if not, the generic C versions are used as a fall-back. On ARM targets,
// the typical path is to have separate builds for NEON and non-NEON devices
// (at least this is true for iOS and Android). Therefore, the NEON support is
// toggled by a build flag: define STBI_NEON to get NEON loops.
//
// The output of the JPEG decoder is slightly different from versions where
// SIMD support was introduced (that is, for versions before 1.49). The
// difference is only +-1 in the 8-bit RGB channels, and only on a small
// fraction of pixels. You can force the pre-1.49 behavior by defining
// STBI_JPEG_OLD, but this will disable some of the SIMD decoding path
// and hence cost some performance.
//
// If for some reason you do not want to use any of SIMD code, or if
// you have issues compiling it, you can disable it entirely by
// defining STBI_NO_SIMD.
//
// ===========================================================================
//
// HDR image support   (disable by defining STBI_NO_HDR)
//
// stb_image now supports loading HDR images in general, and currently
// the Radiance .HDR file format, although the support is provided
// generically. You can still load any file through the existing interface;
// if you attempt to load an HDR file, it will be automatically remapped to
// LDR, assuming gamma 2.2 and an arbitrary scale factor defaulting to 1;
// both of these constants can be reconfigured through this interface:
//
//     stbi_hdr_to_ldr_gamma(2.2f);
//     stbi_hdr_to_ldr_scale(1.0f);
//
// (note, do not use _inverse_ constants; stbi_image will invert them
// appropriately).
//
// Additionally, there is a new, parallel interface for loading files as
// (linear) floats to preserve the full dynamic range:
//
//    float *data = stbi_loadf(filename, &x, &y, &n, 0);
//
// If you load LDR images through this interface, those images will
// be promoted to floating point values, run through the inverse of
// constants corresponding to the above:
//
//     stbi_ldr_to_hdr_scale(1.0f);
//     stbi_ldr_to_hdr_gamma(2.2f);
//
// Finally, given a filename (or an open file or memory block--see header
// file for details) containing image data, you can query for the "most
// appropriate" interface to use (that is, whether the image is HDR or
// not), using:
//
//     stbi_is_hdr(char *filename);
//
// ===========================================================================
//
// iPhone PNG support:
//
// By default we convert iphone-formatted PNGs back to RGB, even though
// they are internally encoded differently. You can disable this conversion
// by by calling stbi_convert_iphone_png_to_rgb(0), in which case
// you will always just get the native iphone "format" through (which
// is BGR stored in RGB).
//
// Call stbi_set_unpremultiply_on_load(1) as well to force a divide per
// pixel to remove any premultiplied alpha *only* if the image file explicitly
// says there's premultiplied data (currently only happens in iPhone images,
// and only if iPhone convert-to-rgb processing is on).
//


#ifndef STBI_NO_STDIO
#include <stdio.h>
#endif // STBI_NO_STDIO

#define STBI_VERSION 1

enum
{
   STBI_default = 0, // only used for req_comp

   STBI_grey       = 1,
   STBI_grey_alpha = 2,
   STBI_rgb        = 3,
   STBI_rgb_alpha  = 4
};

//////////////////////////////////////////////////////////////////////////////
//
// PRIMARY API - works on images of any type
//

//
// load image by filename, open file, or memory buffer
//

typedef struct
{
   int      (*read)  (void *user,char *data,int size);   // fill 'data' with 'size' bytes.  return number of bytes actually read
   void     (*skip)  (void *user,int n);                 // skip the next 'n' bytes, or 'unget' the last -n bytes if negative
   int      (*eof)   (void *user);                       // returns nonzero if we are at end of file/data
} stbi_io_callbacks;

//STBIDEF stbi_uc *stbi_load               (char              const *filename,           int *x, int *y, int *comp, int req_comp);
STBIDEF stbi_uc *stbi_load_from_memory   (stbi_uc           *buffer, int len         , int *x, int *y, int *comp, int req_comp);
STBIDEF stbi_uc *stbi_load_from_callbacks(stbi_io_callbacks *clbk  , void *user      , int *x, int *y, int *comp, int req_comp);

#ifndef STBI_NO_STDIO
STBIDEF stbi_uc *stbi_load_from_file  (FILE *f,                  int *x, int *y, int *comp, int req_comp);
// for stbi_load_from_file, file pointer is left pointing immediately after image
#endif

#ifndef STBI_NO_LINEAR
   STBIDEF float *stbi_loadf                 (char const *filename,           int *x, int *y, int *comp, int req_comp);
   STBIDEF float *stbi_loadf_from_memory     (stbi_uc *buffer, int len, int *x, int *y, int *comp, int req_comp);
   STBIDEF float *stbi_loadf_from_callbacks  (stbi_io_callbacks *clbk, void *user, int *x, int *y, int *comp, int req_comp);

   #ifndef STBI_NO_STDIO
   STBIDEF float *stbi_loadf_from_file  (FILE *f,                int *x, int *y, int *comp, int req_comp);
   #endif
#endif

#ifndef STBI_NO_HDR
   STBIDEF void   stbi_hdr_to_ldr_gamma(float gamma);
   STBIDEF void   stbi_hdr_to_ldr_scale(float scale);
#endif

#ifndef STBI_NO_LINEAR
   STBIDEF void   stbi_ldr_to_hdr_gamma(float gamma);
   STBIDEF void   stbi_ldr_to_hdr_scale(float scale);
#endif // STBI_NO_HDR

// stbi_is_hdr is always defined, but always returns false if STBI_NO_HDR
STBIDEF int    stbi_is_hdr_from_callbacks(stbi_io_callbacks *clbk, void *user);
STBIDEF int    stbi_is_hdr_from_memory(stbi_uc *buffer, int len);
#ifndef STBI_NO_STDIO
STBIDEF int      stbi_is_hdr          (char const *filename);
STBIDEF int      stbi_is_hdr_from_file(FILE *f);
#endif // STBI_NO_STDIO


// get image dimensions & components without fully decoding
STBIDEF int      stbi_info_from_memory(stbi_uc *buffer, int len, int *x, int *y, int *comp);
STBIDEF int      stbi_info_from_callbacks(stbi_io_callbacks *clbk, void *user, int *x, int *y, int *comp);

#ifndef STBI_NO_STDIO
STBIDEF int      stbi_info            (char const *filename,     int *x, int *y, int *comp);
STBIDEF int      stbi_info_from_file  (FILE *f,                  int *x, int *y, int *comp);

#endif



// for image formats that explicitly notate that they have premultiplied alpha,
// we just return the colors as stored in the file. set this flag to force
// unpremultiplication. results are undefined if the unpremultiply overflow.
STBIDEF void stbi_set_unpremultiply_on_load(int flag_true_if_should_unpremultiply);

// indicate whether we should process iphone images back to canonical format,
// or just pass them through "as-is"
STBIDEF void stbi_convert_iphone_png_to_rgb(int flag_true_if_should_convert);

// flip the image vertically, so the first pixel in the output array is the bottom left
STBIDEF void stbi_set_flip_vertically_on_load(int flag_true_if_should_flip);

// ZLIB client - used by PNG, available for other purposes

STBIDEF char *stbi_zlib_decode_malloc_guesssize(const char *buffer, int len, int initial_size, int *outlen);
STBIDEF char *stbi_zlib_decode_malloc_guesssize_headerflag(const char *buffer, int len, int initial_size, int *outlen, int parse_header);
STBIDEF char *stbi_zlib_decode_malloc(const char *buffer, int len, int *outlen);
STBIDEF int   stbi_zlib_decode_buffer(char *obuffer, int olen, const char *ibuffer, int ilen);

STBIDEF char *stbi_zlib_decode_noheader_malloc(const char *buffer, int len, int *outlen);
STBIDEF int   stbi_zlib_decode_noheader_buffer(char *obuffer, int olen, const char *ibuffer, int ilen);


#if defined(STBI_ONLY_JPEG) || defined(STBI_ONLY_PNG) || defined(STBI_ONLY_BMP) \
  || defined(STBI_ONLY_TGA) || defined(STBI_ONLY_GIF) || defined(STBI_ONLY_PSD) \
  || defined(STBI_ONLY_HDR) || defined(STBI_ONLY_PIC) || defined(STBI_ONLY_PNM) \
  || defined(STBI_ONLY_ZLIB)
   #ifndef STBI_ONLY_JPEG
   #define STBI_NO_JPEG
   #endif
   #ifndef STBI_ONLY_PNG
   #define STBI_NO_PNG
   #endif
   #ifndef STBI_ONLY_BMP
   #define STBI_NO_BMP
   #endif
   #ifndef STBI_ONLY_PSD
   #define STBI_NO_PSD
   #endif
   #ifndef STBI_ONLY_TGA
   #define STBI_NO_TGA
   #endif
   #ifndef STBI_ONLY_GIF
   #define STBI_NO_GIF
   #endif
   #ifndef STBI_ONLY_HDR
   #define STBI_NO_HDR
   #endif
   #ifndef STBI_ONLY_PIC
   #define STBI_NO_PIC
   #endif
   #ifndef STBI_ONLY_PNM
   #define STBI_NO_PNM
   #endif
#endif

#if defined(STBI_NO_PNG) && !defined(STBI_SUPPORT_ZLIB) && !defined(STBI_NO_ZLIB)
#define STBI_NO_ZLIB
#endif


#include <stdarg.h>
#include <stddef.h> // ptrdiff_t on osx
#include <stdlib.h>
#include <string.h>

#if !defined(STBI_NO_LINEAR) || !defined(STBI_NO_HDR)
#include <math.h>  // ldexp
#endif

#ifndef STBI_NO_STDIO
#include <stdio.h>
#endif

#ifndef STBI_ASSERT
#include <assert.h>
#define STBI_ASSERT(x) assert(x)
#endif


#ifndef _MSC_VER
   #define stbi_inline inline
#else
   #define stbi_inline __forceinline
#endif


#ifdef _MSC_VER
typedef unsigned short stbi__uint16;
typedef   signed short stbi__int16;
typedef unsigned int   stbi__uint32;
typedef   signed int   stbi__int32;
#else
#include <stdint.h>
typedef uint16_t stbi__uint16;
typedef int16_t  stbi__int16;
typedef uint32_t stbi__uint32;
typedef int32_t  stbi__int32;
#endif

// should produce compiler error if size is wrong
typedef unsigned char validate_uint32[sizeof(stbi__uint32)==4 ? 1 : -1];

#ifdef _MSC_VER
#define STBI_NOTUSED(v)  (void)(v)
#else
#define STBI_NOTUSED(v)  (void)sizeof(v)
#endif

#ifdef _MSC_VER
#define STBI_HAS_LROTL
#endif

#ifdef STBI_HAS_LROTL
   #define stbi_lrot(x,y)  _lrotl(x,y)
#else
   #define stbi_lrot(x,y)  (((x) << (y)) | ((x) >> (32 - (y))))
#endif

#if defined(STBI_MALLOC) && defined(STBI_FREE) && defined(STBI_REALLOC)
// ok
#elif !defined(STBI_MALLOC) && !defined(STBI_FREE) && !defined(STBI_REALLOC)
// ok
#else
#error "Must define all or none of STBI_MALLOC, STBI_FREE, and STBI_REALLOC."
#endif

#ifndef STBI_MALLOC
#define STBI_MALLOC(sz)    malloc(sz)
#define STBI_REALLOC(p,sz) realloc(p,sz)
#define STBI_FREE(p)       free(p)
#endif

// x86/x64 detection
#if defined(__x86_64__) || defined(_M_X64)
#define STBI__X64_TARGET
#elif defined(__i386) || defined(_M_IX86)
#define STBI__X86_TARGET
#endif

#if defined(__GNUC__) && (defined(STBI__X86_TARGET) || defined(STBI__X64_TARGET)) && !defined(__SSE2__) && !defined(STBI_NO_SIMD)
// NOTE: not clear do we actually need this for the 64-bit path?
// gcc doesn't support sse2 intrinsics unless you compile with -msse2,
// (but compiling with -msse2 allows the compiler to use SSE2 everywhere;
// this is just broken and gcc are jerks for not fixing it properly
// http://www.virtualdub.org/blog/pivot/entry.php?id=363 )
#define STBI_NO_SIMD
#endif

#if defined(__MINGW32__) && defined(STBI__X86_TARGET) && !defined(STBI_MINGW_ENABLE_SSE2) && !defined(STBI_NO_SIMD)
// Note that __MINGW32__ doesn't actually mean 32-bit, so we have to avoid STBI__X64_TARGET
//
// 32-bit MinGW wants ESP to be 16-byte aligned, but this is not in the
// Windows ABI and VC++ as well as Windows DLLs don't maintain that invariant.
// As a result, enabling SSE2 on 32-bit MinGW is dangerous when not
// simultaneously enabling "-mstackrealign".
//
// See https://github.com/nothings/stb/issues/81 for more information.
//
// So default to no SSE2 on 32-bit MinGW. If you've read this far and added
// -mstackrealign to your build settings, feel free to #define STBI_MINGW_ENABLE_SSE2.
#define STBI_NO_SIMD
#endif

#if !defined(STBI_NO_SIMD) && defined(STBI__X86_TARGET)
#define STBI_SSE2
#include <emmintrin.h>

#ifdef _MSC_VER

#if _MSC_VER >= 1400  // not VC6
#include <intrin.h> // __cpuid
static int stbi__cpuid3(void)
{
   int info[4];
   __cpuid(info,1);
   return info[3];
}
#else
static int stbi__cpuid3(void)
{
   int res;
   __asm {
      mov  eax,1
      cpuid
      mov  res,edx
   }
   return res;
}
#endif

#define STBI_SIMD_ALIGN(type, name) __declspec(align(16)) type name

static int stbi__sse2_available()
{
   int info3 = stbi__cpuid3();
   return ((info3 >> 26) & 1) != 0;
}
#else // assume GCC-style if not VC++
#define STBI_SIMD_ALIGN(type, name) type name __attribute__((aligned(16)))

static int stbi__sse2_available()
{
#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 // GCC 4.8 or later
   // GCC 4.8+ has a nice way to do this
   return __builtin_cpu_supports("sse2");
#else
   // portable way to do this, preferably without using GCC inline ASM?
   // just bail for now.
   return 0;
#endif
}
#endif
#endif

// ARM NEON
#if defined(STBI_NO_SIMD) && defined(STBI_NEON)
#undef STBI_NEON
#endif

#ifdef STBI_NEON
#include <arm_neon.h>
// assume GCC or Clang on ARM targets
#define STBI_SIMD_ALIGN(type, name) type name __attribute__((aligned(16)))
#endif

#ifndef STBI_SIMD_ALIGN
#define STBI_SIMD_ALIGN(type, name) type name
#endif

///////////////////////////////////////////////
//
//  stbi__context struct and start_xxx functions

// stbi__context structure is our basic context used by all images, so it
// contains all the IO context, plus some basic image information
typedef struct
{
   stbi__uint32 img_x, img_y;
   int img_n, img_out_n;

   stbi_io_callbacks io;
   void *io_user_data;

   int read_from_callbacks;
   int buflen;
   stbi_uc buffer_start[128];

   stbi_uc *img_buffer, *img_buffer_end;
   stbi_uc *img_buffer_original;
} stbi__context;


static void stbi__refill_buffer(stbi__context *s);

// initialize a memory-decode context
static void stbi__start_mem(stbi__context *s, stbi_uc *buffer, int len)
{
   s->io.read = NULL;
   s->read_from_callbacks = 0;
   s->img_buffer = s->img_buffer_original = buffer;
   s->img_buffer_end = buffer+len;
}

// initialize a callback-based context
static void stbi__start_callbacks(stbi__context *s, stbi_io_callbacks *c, void *user)
{
   s->io = *c;
   s->io_user_data = user;
   s->buflen = sizeof(s->buffer_start);
   s->read_from_callbacks = 1;
   s->img_buffer_original = s->buffer_start;
   stbi__refill_buffer(s);
}

#ifndef STBI_NO_STDIO

static int stbi__stdio_read(void *user, char *data, int size)
{
   return (int) fread(data,1,size,(FILE*) user);
}

static void stbi__stdio_skip(void *user, int n)
{
   fseek((FILE*) user, n, SEEK_CUR);
}

static int stbi__stdio_eof(void *user)
{
   return feof((FILE*) user);
}

static stbi_io_callbacks stbi__stdio_callbacks =
{
   stbi__stdio_read,
   stbi__stdio_skip,
   stbi__stdio_eof,
};

static void stbi__start_file(stbi__context *s, FILE *f)
{
   stbi__start_callbacks(s, &stbi__stdio_callbacks, (void *) f);
}

//static void stop_file(stbi__context *s) { }

#endif // !STBI_NO_STDIO

static void stbi__rewind(stbi__context *s)
{
   // conceptually rewind SHOULD rewind to the beginning of the stream,
   // but we just rewind to the beginning of the initial buffer, because
   // we only use it after doing 'test', which only ever looks at at most 92 bytes
   s->img_buffer = s->img_buffer_original;
}

#ifndef STBI_NO_JPEG
static int      stbi__jpeg_test(stbi__context *s);
static stbi_uc *stbi__jpeg_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__jpeg_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_PNG
static int      stbi__png_test(stbi__context *s);
static stbi_uc *stbi__png_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__png_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_BMP
static int      stbi__bmp_test(stbi__context *s);
static stbi_uc *stbi__bmp_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__bmp_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_TGA
static int      stbi__tga_test(stbi__context *s);
static stbi_uc *stbi__tga_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__tga_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_PSD
static int      stbi__psd_test(stbi__context *s);
static stbi_uc *stbi__psd_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__psd_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_HDR
static int      stbi__hdr_test(stbi__context *s);
static float   *stbi__hdr_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__hdr_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_PIC
static int      stbi__pic_test(stbi__context *s);
static stbi_uc *stbi__pic_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__pic_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_GIF
static int      stbi__gif_test(stbi__context *s);
static stbi_uc *stbi__gif_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__gif_info(stbi__context *s, int *x, int *y, int *comp);
#endif

#ifndef STBI_NO_PNM
static int      stbi__pnm_test(stbi__context *s);
static stbi_uc *stbi__pnm_load(stbi__context *s, int *x, int *y, int *comp, int req_comp);
static int      stbi__pnm_info(stbi__context *s, int *x, int *y, int *comp);
#endif

// this is not threadsafe
static const char *stbi__g_failure_reason;

STBIDEF const char *stbi_failure_reason(void)
{
   return stbi__g_failure_reason;
}

static int stbi__err(const char *str)
{
   stbi__g_failure_reason = str;
   return 0;
}

static void *stbi__malloc(size_t size)
{
    return STBI_MALLOC(size);
}

// stbi__err - error
// stbi__errpf - error returning pointer to float
// stbi__errpuc - error returning pointer to unsigned char

#ifdef STBI_NO_FAILURE_STRINGS
   #define stbi__err(x,y)  0
#elif defined(STBI_FAILURE_USERMSG)
   #define stbi__err(x,y)  stbi__err(y)
#else
   #define stbi__err(x,y)  stbi__err(x)
#endif

#define stbi__errpf(x,y)   ((float *) (stbi__err(x,y)?NULL:NULL))
#define stbi__errpuc(x,y)  ((unsigned char *) (stbi__err(x,y)?NULL:NULL))

STBIDEF void stbi_image_free(void *retval_from_stbi_load)
{
   STBI_FREE(retval_from_stbi_load);
}

#ifndef STBI_NO_LINEAR
static float   *stbi__ldr_to_hdr(stbi_uc *data, int x, int y, int comp);
#endif

#ifndef STBI_NO_HDR
static stbi_uc *stbi__hdr_to_ldr(float   *data, int x, int y, int comp);
#endif

static int stbi__vertically_flip_on_load = 0;

STBIDEF void stbi_set_flip_vertically_on_load(int flag_true_if_should_flip)
{
    stbi__vertically_flip_on_load = flag_true_if_should_flip;
}

static unsigned char *stbi__load_main(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   #ifndef STBI_NO_JPEG
   if (stbi__jpeg_test(s)) return stbi__jpeg_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_PNG
   if (stbi__png_test(s))  return stbi__png_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_BMP
   if (stbi__bmp_test(s))  return stbi__bmp_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_GIF
   if (stbi__gif_test(s))  return stbi__gif_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_PSD
   if (stbi__psd_test(s))  return stbi__psd_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_PIC
   if (stbi__pic_test(s))  return stbi__pic_load(s,x,y,comp,req_comp);
   #endif
   #ifndef STBI_NO_PNM
   if (stbi__pnm_test(s))  return stbi__pnm_load(s,x,y,comp,req_comp);
   #endif

   #ifndef STBI_NO_HDR
   if (stbi__hdr_test(s)) {
      float *hdr = stbi__hdr_load(s, x,y,comp,req_comp);
      return stbi__hdr_to_ldr(hdr, *x, *y, req_comp ? req_comp : *comp);
   }
   #endif

   #ifndef STBI_NO_TGA
   // test tga last because it's a crappy test!
   if (stbi__tga_test(s))
      return stbi__tga_load(s,x,y,comp,req_comp);
   #endif

   return stbi__errpuc("unknown image type", "Image not of any known type, or corrupt");
}

static unsigned char *stbi__load_flip(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   unsigned char *result = stbi__load_main(s, x, y, comp, req_comp);

   if (stbi__vertically_flip_on_load && result != NULL) {
      int w = *x, h = *y;
      int depth = req_comp ? req_comp : *comp;
      int row,col,z;
      stbi_uc temp;

      // @OPTIMIZE: use a bigger temp buffer and memcpy multiple pixels at once
      for (row = 0; row < (h>>1); row++) {
         for (col = 0; col < w; col++) {
            for (z = 0; z < depth; z++) {
               temp = result[(row * w + col) * depth + z];
               result[(row * w + col) * depth + z] = result[((h - row - 1) * w + col) * depth + z];
               result[((h - row - 1) * w + col) * depth + z] = temp;
            }
         }
      }
   }

   return result;
}

static void stbi__float_postprocess(float *result, int *x, int *y, int *comp, int req_comp)
{
   if (stbi__vertically_flip_on_load && result != NULL) {
      int w = *x, h = *y;
      int depth = req_comp ? req_comp : *comp;
      int row,col,z;
      float temp;

      // @OPTIMIZE: use a bigger temp buffer and memcpy multiple pixels at once
      for (row = 0; row < (h>>1); row++) {
         for (col = 0; col < w; col++) {
            for (z = 0; z < depth; z++) {
               temp = result[(row * w + col) * depth + z];
               result[(row * w + col) * depth + z] = result[((h - row - 1) * w + col) * depth + z];
               result[((h - row - 1) * w + col) * depth + z] = temp;
            }
         }
      }
   }
}


#ifndef STBI_NO_STDIO

static FILE *stbi__fopen(char const *filename, char const *mode)
{
   FILE *f;
#if defined(_MSC_VER) && _MSC_VER >= 1400
   if (0 != fopen_s(&f, filename, mode))
      f=0;
#else
   f = fopen(filename, mode);
#endif
   return f;
}


stbi_uc *stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp)
{
   FILE *f = stbi__fopen(filename, "rb");
   unsigned char *result;
   if (!f) return stbi__errpuc("can't fopen", "Unable to open file");
   result = stbi_load_from_file(f,x,y,comp,req_comp);
   fclose(f);
   return result;
}

STBIDEF stbi_uc *stbi_load_from_file(FILE *f, int *x, int *y, int *comp, int req_comp)
{
   unsigned char *result;
   stbi__context s;
   stbi__start_file(&s,f);
   result = stbi__load_flip(&s,x,y,comp,req_comp);
   if (result) {
      // need to 'unget' all the characters in the IO buffer
      fseek(f, - (int) (s.img_buffer_end - s.img_buffer), SEEK_CUR);
   }
   return result;
}
#endif //!STBI_NO_STDIO

STBIDEF stbi_uc *stbi_load_from_memory(stbi_uc *buffer, int len, int *x, int *y, int *comp, int req_comp)
{
   stbi__context s;
   stbi__start_mem(&s,buffer,len);
   return stbi__load_flip(&s,x,y,comp,req_comp);
}

STBIDEF stbi_uc *stbi_load_from_callbacks(stbi_io_callbacks *clbk, void *user, int *x, int *y, int *comp, int req_comp)
{
   stbi__context s;
   stbi__start_callbacks(&s, clbk, user);
   return stbi__load_flip(&s,x,y,comp,req_comp);
}

#ifndef STBI_NO_LINEAR
static float *stbi__loadf_main(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   unsigned char *data;
   #ifndef STBI_NO_HDR
   if (stbi__hdr_test(s)) {
      float *hdr_data = stbi__hdr_load(s,x,y,comp,req_comp);
      if (hdr_data)
         stbi__float_postprocess(hdr_data,x,y,comp,req_comp);
      return hdr_data;
   }
   #endif
   data = stbi__load_flip(s, x, y, comp, req_comp);
   if (data)
      return stbi__ldr_to_hdr(data, *x, *y, req_comp ? req_comp : *comp);
   return stbi__errpf("unknown image type", "Image not of any known type, or corrupt");
}

STBIDEF float *stbi_loadf_from_memory(stbi_uc *buffer, int len, int *x, int *y, int *comp, int req_comp)
{
   stbi__context s;
   stbi__start_mem(&s,buffer,len);
   return stbi__loadf_main(&s,x,y,comp,req_comp);
}

STBIDEF float *stbi_loadf_from_callbacks(stbi_io_callbacks *clbk, void *user, int *x, int *y, int *comp, int req_comp)
{
   stbi__context s;
   stbi__start_callbacks(&s, clbk, user);
   return stbi__loadf_main(&s,x,y,comp,req_comp);
}

#ifndef STBI_NO_STDIO
STBIDEF float *stbi_loadf(char const *filename, int *x, int *y, int *comp, int req_comp)
{
   float *result;
   FILE *f = stbi__fopen(filename, "rb");
   if (!f) return stbi__errpf("can't fopen", "Unable to open file");
   result = stbi_loadf_from_file(f,x,y,comp,req_comp);
   fclose(f);
   return result;
}

STBIDEF float *stbi_loadf_from_file(FILE *f, int *x, int *y, int *comp, int req_comp)
{
   stbi__context s;
   stbi__start_file(&s,f);
   return stbi__loadf_main(&s,x,y,comp,req_comp);
}
#endif // !STBI_NO_STDIO

#endif // !STBI_NO_LINEAR

// these is-hdr-or-not is defined independent of whether STBI_NO_LINEAR is
// defined, for API simplicity; if STBI_NO_LINEAR is defined, it always
// reports false!

STBIDEF int stbi_is_hdr_from_memory(stbi_uc *buffer, int len)
{
   #ifndef STBI_NO_HDR
   stbi__context s;
   stbi__start_mem(&s,buffer,len);
   return stbi__hdr_test(&s);
   #else
   STBI_NOTUSED(buffer);
   STBI_NOTUSED(len);
   return 0;
   #endif
}

#ifndef STBI_NO_STDIO
STBIDEF int      stbi_is_hdr          (char const *filename)
{
   FILE *f = stbi__fopen(filename, "rb");
   int result=0;
   if (f) {
      result = stbi_is_hdr_from_file(f);
      fclose(f);
   }
   return result;
}

STBIDEF int      stbi_is_hdr_from_file(FILE *f)
{
   #ifndef STBI_NO_HDR
   stbi__context s;
   stbi__start_file(&s,f);
   return stbi__hdr_test(&s);
   #else
   return 0;
   #endif
}
#endif // !STBI_NO_STDIO

STBIDEF int      stbi_is_hdr_from_callbacks(stbi_io_callbacks *clbk, void *user)
{
   #ifndef STBI_NO_HDR
   stbi__context s;
   stbi__start_callbacks(&s, clbk, user);
   return stbi__hdr_test(&s);
   #else
   return 0;
   #endif
}

static float stbi__h2l_gamma_i=1.0f/2.2f, stbi__h2l_scale_i=1.0f;
static float stbi__l2h_gamma=2.2f, stbi__l2h_scale=1.0f;

#ifndef STBI_NO_LINEAR
STBIDEF void   stbi_ldr_to_hdr_gamma(float gamma) { stbi__l2h_gamma = gamma; }
STBIDEF void   stbi_ldr_to_hdr_scale(float scale) { stbi__l2h_scale = scale; }
#endif

STBIDEF void   stbi_hdr_to_ldr_gamma(float gamma) { stbi__h2l_gamma_i = 1/gamma; }
STBIDEF void   stbi_hdr_to_ldr_scale(float scale) { stbi__h2l_scale_i = 1/scale; }


//////////////////////////////////////////////////////////////////////////////
//
// Common code used by all image loaders
//

enum
{
   STBI__SCAN_load=0,
   STBI__SCAN_type,
   STBI__SCAN_header
};

static void stbi__refill_buffer(stbi__context *s)
{
   int n = (s->io.read)(s->io_user_data,(char*)s->buffer_start,s->buflen);
   if (n == 0) {
      // at end of file, treat same as if from memory, but need to handle case
      // where s->img_buffer isn't pointing to safe memory, e.g. 0-byte file
      s->read_from_callbacks = 0;
      s->img_buffer = s->buffer_start;
      s->img_buffer_end = s->buffer_start+1;
      *s->img_buffer = 0;
   } else {
      s->img_buffer = s->buffer_start;
      s->img_buffer_end = s->buffer_start + n;
   }
}

stbi_inline static stbi_uc stbi__get8(stbi__context *s)
{
   if (s->img_buffer < s->img_buffer_end)
      return *s->img_buffer++;
   if (s->read_from_callbacks) {
      stbi__refill_buffer(s);
      return *s->img_buffer++;
   }
   return 0;
}

stbi_inline static int stbi__at_eof(stbi__context *s)
{
   if (s->io.read) {
      if (!(s->io.eof)(s->io_user_data)) return 0;
      // if feof() is true, check if buffer = end
      // special case: we've only got the special 0 character at the end
      if (s->read_from_callbacks == 0) return 1;
   }

   return s->img_buffer >= s->img_buffer_end;
}

static void stbi__skip(stbi__context *s, int n)
{
   if (n < 0) {
      s->img_buffer = s->img_buffer_end;
      return;
   }
   if (s->io.read) {
      int blen = (int) (s->img_buffer_end - s->img_buffer);
      if (blen < n) {
         s->img_buffer = s->img_buffer_end;
         (s->io.skip)(s->io_user_data, n - blen);
         return;
      }
   }
   s->img_buffer += n;
}

static int stbi__getn(stbi__context *s, stbi_uc *buffer, int n)
{
   if (s->io.read) {
      int blen = (int) (s->img_buffer_end - s->img_buffer);
      if (blen < n) {
         int res, count;

         memcpy(buffer, s->img_buffer, blen);

         count = (s->io.read)(s->io_user_data, (char*) buffer + blen, n - blen);
         res = (count == (n-blen));
         s->img_buffer = s->img_buffer_end;
         return res;
      }
   }

   if (s->img_buffer+n <= s->img_buffer_end) {
      memcpy(buffer, s->img_buffer, n);
      s->img_buffer += n;
      return 1;
   } else
      return 0;
}

static int stbi__get16be(stbi__context *s)
{
   int z = stbi__get8(s);
   return (z << 8) + stbi__get8(s);
}

static stbi__uint32 stbi__get32be(stbi__context *s)
{
   stbi__uint32 z = stbi__get16be(s);
   return (z << 16) + stbi__get16be(s);
}

static int stbi__get16le(stbi__context *s)
{
   int z = stbi__get8(s);
   return z + (stbi__get8(s) << 8);
}

static stbi__uint32 stbi__get32le(stbi__context *s)
{
   stbi__uint32 z = stbi__get16le(s);
   return z + (stbi__get16le(s) << 16);
}

#define STBI__BYTECAST(x)  ((stbi_uc) ((x) & 255))  // truncate int to byte without warnings


//////////////////////////////////////////////////////////////////////////////
//
//  generic converter from built-in img_n to req_comp
//    individual types do this automatically as much as possible (e.g. jpeg
//    does all cases internally since it needs to colorspace convert anyway,
//    and it never has alpha, so very few cases ). png can automatically
//    interleave an alpha=255 channel, but falls back to this for other cases
//
//  assume data buffer is malloced, so malloc a new one and free that one
//  only failure mode is malloc failing

static stbi_uc stbi__compute_y(int r, int g, int b)
{
   return (stbi_uc) (((r*77) + (g*150) +  (29*b)) >> 8);
}

static unsigned char *stbi__convert_format(unsigned char *data, int img_n, int req_comp, unsigned int x, unsigned int y)
{
   int i,j;
   unsigned char *good;

   if (req_comp == img_n) return data;
   STBI_ASSERT(req_comp >= 1 && req_comp <= 4);

   good = (unsigned char *) stbi__malloc(req_comp * x * y);
   if (good == NULL) {
      STBI_FREE(data);
      return stbi__errpuc("outofmem", "Out of memory");
   }

   for (j=0; j < (int) y; ++j) {
      unsigned char *src  = data + j * x * img_n   ;
      unsigned char *dest = good + j * x * req_comp;

      #define COMBO(a,b)  ((a)*8+(b))
      #define CASE(a,b)   case COMBO(a,b): for(i=x-1; i >= 0; --i, src += a, dest += b)
      // convert source image with img_n components to one with req_comp components;
      // avoid switch per pixel, so use switch per scanline and massive macros
      switch (COMBO(img_n, req_comp)) {
         CASE(1,2) dest[0]=src[0], dest[1]=255; break;
         CASE(1,3) dest[0]=dest[1]=dest[2]=src[0]; break;
         CASE(1,4) dest[0]=dest[1]=dest[2]=src[0], dest[3]=255; break;
         CASE(2,1) dest[0]=src[0]; break;
         CASE(2,3) dest[0]=dest[1]=dest[2]=src[0]; break;
         CASE(2,4) dest[0]=dest[1]=dest[2]=src[0], dest[3]=src[1]; break;
         CASE(3,4) dest[0]=src[0],dest[1]=src[1],dest[2]=src[2],dest[3]=255; break;
         CASE(3,1) dest[0]=stbi__compute_y(src[0],src[1],src[2]); break;
         CASE(3,2) dest[0]=stbi__compute_y(src[0],src[1],src[2]), dest[1] = 255; break;
         CASE(4,1) dest[0]=stbi__compute_y(src[0],src[1],src[2]); break;
         CASE(4,2) dest[0]=stbi__compute_y(src[0],src[1],src[2]), dest[1] = src[3]; break;
         CASE(4,3) dest[0]=src[0],dest[1]=src[1],dest[2]=src[2]; break;
         default: STBI_ASSERT(0);
      }
      #undef CASE
   }

   STBI_FREE(data);
   return good;
}

#ifndef STBI_NO_LINEAR
static float   *stbi__ldr_to_hdr(stbi_uc *data, int x, int y, int comp)
{
   int i,k,n;
   float *output = (float *) stbi__malloc(x * y * comp * sizeof(float));
   if (output == NULL) { STBI_FREE(data); return stbi__errpf("outofmem", "Out of memory"); }
   // compute number of non-alpha components
   if (comp & 1) n = comp; else n = comp-1;
   for (i=0; i < x*y; ++i) {
      for (k=0; k < n; ++k) {
         output[i*comp + k] = (float) (pow(data[i*comp+k]/255.0f, stbi__l2h_gamma) * stbi__l2h_scale);
      }
      if (k < comp) output[i*comp + k] = data[i*comp+k]/255.0f;
   }
   STBI_FREE(data);
   return output;
}
#endif

#ifndef STBI_NO_HDR
#define stbi__float2int(x)   ((int) (x))
static stbi_uc *stbi__hdr_to_ldr(float   *data, int x, int y, int comp)
{
   int i,k,n;
   stbi_uc *output = (stbi_uc *) stbi__malloc(x * y * comp);
   if (output == NULL) { STBI_FREE(data); return stbi__errpuc("outofmem", "Out of memory"); }
   // compute number of non-alpha components
   if (comp & 1) n = comp; else n = comp-1;
   for (i=0; i < x*y; ++i) {
      for (k=0; k < n; ++k) {
         float z = (float) pow(data[i*comp+k]*stbi__h2l_scale_i, stbi__h2l_gamma_i) * 255 + 0.5f;
         if (z < 0) z = 0;
         if (z > 255) z = 255;
         output[i*comp + k] = (stbi_uc) stbi__float2int(z);
      }
      if (k < comp) {
         float z = data[i*comp+k] * 255 + 0.5f;
         if (z < 0) z = 0;
         if (z > 255) z = 255;
         output[i*comp + k] = (stbi_uc) stbi__float2int(z);
      }
   }
   STBI_FREE(data);
   return output;
}
#endif

//////////////////////////////////////////////////////////////////////////////
//
//  "baseline" JPEG/JFIF decoder
//
//    simple implementation
//      - doesn't support delayed output of y-dimension
//      - simple interface (only one output format: 8-bit interleaved RGB)
//      - doesn't try to recover corrupt jpegs
//      - doesn't allow partial loading, loading multiple at once
//      - still fast on x86 (copying globals into locals doesn't help x86)
//      - allocates lots of intermediate memory (full size of all components)
//        - non-interleaved case requires this anyway
//        - allows good upsampling (see next)
//    high-quality
//      - upsampled channels are bilinearly interpolated, even across blocks
//      - quality integer IDCT derived from IJG's 'slow'
//    performance
//      - fast huffman; reasonable integer IDCT
//      - some SIMD kernels for common paths on targets with SSE2/NEON
//      - uses a lot of intermediate memory, could cache poorly

#ifndef STBI_NO_JPEG

// huffman decoding acceleration
#define FAST_BITS   9  // larger handles more cases; smaller stomps less cache

typedef struct
{
   stbi_uc  fast[1 << FAST_BITS];
   // weirdly, repacking this into AoS is a 10% speed loss, instead of a win
   stbi__uint16 code[256];
   stbi_uc  values[256];
   stbi_uc  size[257];
   unsigned int maxcode[18];
   int    delta[17];   // old 'firstsymbol' - old 'firstcode'
} stbi__huffman;

typedef struct
{
   stbi__context *s;
   stbi__huffman huff_dc[4];
   stbi__huffman huff_ac[4];
   stbi_uc dequant[4][64];
   stbi__int16 fast_ac[4][1 << FAST_BITS];

// sizes for components, interleaved MCUs
   int img_h_max, img_v_max;
   int img_mcu_x, img_mcu_y;
   int img_mcu_w, img_mcu_h;

// definition of jpeg image component
   struct
   {
      int id;
      int h,v;
      int tq;
      int hd,ha;
      int dc_pred;

      int x,y,w2,h2;
      stbi_uc *data;
      void *raw_data, *raw_coeff;
      stbi_uc *linebuf;
      short   *coeff;   // progressive only
      int      coeff_w, coeff_h; // number of 8x8 coefficient blocks
   } img_comp[4];

   stbi__uint32   code_buffer; // jpeg entropy-coded buffer
   int            code_bits;   // number of valid bits
   unsigned char  marker;      // marker seen while filling entropy buffer
   int            nomore;      // flag if we saw a marker so must stop

   int            progressive;
   int            spec_start;
   int            spec_end;
   int            succ_high;
   int            succ_low;
   int            eob_run;

   int scan_n, order[4];
   int restart_interval, todo;

// kernels
   void (*idct_block_kernel)(stbi_uc *out, int out_stride, short data[64]);
   void (*YCbCr_to_RGB_kernel)(stbi_uc *out, const stbi_uc *y, const stbi_uc *pcb, const stbi_uc *pcr, int count, int step);
   stbi_uc *(*resample_row_hv_2_kernel)(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs);
} stbi__jpeg;

static int stbi__build_huffman(stbi__huffman *h, int *count)
{
   int i,j,k=0,code;
   // build size list for each symbol (from JPEG spec)
   for (i=0; i < 16; ++i)
      for (j=0; j < count[i]; ++j)
         h->size[k++] = (stbi_uc) (i+1);
   h->size[k] = 0;

   // compute actual symbols (from jpeg spec)
   code = 0;
   k = 0;
   for(j=1; j <= 16; ++j) {
      // compute delta to add to code to compute symbol id
      h->delta[j] = k - code;
      if (h->size[k] == j) {
         while (h->size[k] == j)
            h->code[k++] = (stbi__uint16) (code++);
         if (code-1 >= (1 << j)) return stbi__err("bad code lengths","Corrupt JPEG");
      }
      // compute largest code + 1 for this size, preshifted as needed later
      h->maxcode[j] = code << (16-j);
      code <<= 1;
   }
   h->maxcode[j] = 0xffffffff;

   // build non-spec acceleration table; 255 is flag for not-accelerated
   memset(h->fast, 255, 1 << FAST_BITS);
   for (i=0; i < k; ++i) {
      int s = h->size[i];
      if (s <= FAST_BITS) {
         int c = h->code[i] << (FAST_BITS-s);
         int m = 1 << (FAST_BITS-s);
         for (j=0; j < m; ++j) {
            h->fast[c+j] = (stbi_uc) i;
         }
      }
   }
   return 1;
}

// build a table that decodes both magnitude and value of small ACs in
// one go.
static void stbi__build_fast_ac(stbi__int16 *fast_ac, stbi__huffman *h)
{
   int i;
   for (i=0; i < (1 << FAST_BITS); ++i) {
      stbi_uc fast = h->fast[i];
      fast_ac[i] = 0;
      if (fast < 255) {
         int rs = h->values[fast];
         int run = (rs >> 4) & 15;
         int magbits = rs & 15;
         int len = h->size[fast];

         if (magbits && len + magbits <= FAST_BITS) {
            // magnitude code followed by receive_extend code
            int k = ((i << len) & ((1 << FAST_BITS) - 1)) >> (FAST_BITS - magbits);
            int m = 1 << (magbits - 1);
            if (k < m) k += (-1 << magbits) + 1;
            // if the result is small enough, we can fit it in fast_ac table
            if (k >= -128 && k <= 127)
               fast_ac[i] = (stbi__int16) ((k << 8) + (run << 4) + (len + magbits));
         }
      }
   }
}

static void stbi__grow_buffer_unsafe(stbi__jpeg *j)
{
   do {
      int b = j->nomore ? 0 : stbi__get8(j->s);
      if (b == 0xff) {
         int c = stbi__get8(j->s);
         if (c != 0) {
            j->marker = (unsigned char) c;
            j->nomore = 1;
            return;
         }
      }
      j->code_buffer |= b << (24 - j->code_bits);
      j->code_bits += 8;
   } while (j->code_bits <= 24);
}

// (1 << n) - 1
static stbi__uint32 stbi__bmask[17]={0,1,3,7,15,31,63,127,255,511,1023,2047,4095,8191,16383,32767,65535};

// decode a jpeg huffman value from the bitstream
stbi_inline static int stbi__jpeg_huff_decode(stbi__jpeg *j, stbi__huffman *h)
{
   unsigned int temp;
   int c,k;

   if (j->code_bits < 16) stbi__grow_buffer_unsafe(j);

   // look at the top FAST_BITS and determine what symbol ID it is,
   // if the code is <= FAST_BITS
   c = (j->code_buffer >> (32 - FAST_BITS)) & ((1 << FAST_BITS)-1);
   k = h->fast[c];
   if (k < 255) {
      int s = h->size[k];
      if (s > j->code_bits)
         return -1;
      j->code_buffer <<= s;
      j->code_bits -= s;
      return h->values[k];
   }

   // naive test is to shift the code_buffer down so k bits are
   // valid, then test against maxcode. To speed this up, we've
   // preshifted maxcode left so that it has (16-k) 0s at the
   // end; in other words, regardless of the number of bits, it
   // wants to be compared against something shifted to have 16;
   // that way we don't need to shift inside the loop.
   temp = j->code_buffer >> 16;
   for (k=FAST_BITS+1 ; ; ++k)
      if (temp < h->maxcode[k])
         break;
   if (k == 17) {
      // error! code not found
      j->code_bits -= 16;
      return -1;
   }

   if (k > j->code_bits)
      return -1;

   // convert the huffman code to the symbol id
   c = ((j->code_buffer >> (32 - k)) & stbi__bmask[k]) + h->delta[k];
   STBI_ASSERT((((j->code_buffer) >> (32 - h->size[c])) & stbi__bmask[h->size[c]]) == h->code[c]);

   // convert the id to a symbol
   j->code_bits -= k;
   j->code_buffer <<= k;
   return h->values[c];
}

// bias[n] = (-1<<n) + 1
static int const stbi__jbias[16] = {0,-1,-3,-7,-15,-31,-63,-127,-255,-511,-1023,-2047,-4095,-8191,-16383,-32767};

// combined JPEG 'receive' and JPEG 'extend', since baseline
// always extends everything it receives.
stbi_inline static int stbi__extend_receive(stbi__jpeg *j, int n)
{
   unsigned int k;
   int sgn;
   if (j->code_bits < n) stbi__grow_buffer_unsafe(j);

   sgn = (stbi__int32)j->code_buffer >> 31; // sign bit is always in MSB
   k = stbi_lrot(j->code_buffer, n);
   STBI_ASSERT(n >= 0 && n < (int) (sizeof(stbi__bmask)/sizeof(*stbi__bmask)));
   j->code_buffer = k & ~stbi__bmask[n];
   k &= stbi__bmask[n];
   j->code_bits -= n;
   return k + (stbi__jbias[n] & ~sgn);
}

// get some unsigned bits
stbi_inline static int stbi__jpeg_get_bits(stbi__jpeg *j, int n)
{
   unsigned int k;
   if (j->code_bits < n) stbi__grow_buffer_unsafe(j);
   k = stbi_lrot(j->code_buffer, n);
   j->code_buffer = k & ~stbi__bmask[n];
   k &= stbi__bmask[n];
   j->code_bits -= n;
   return k;
}

stbi_inline static int stbi__jpeg_get_bit(stbi__jpeg *j)
{
   unsigned int k;
   if (j->code_bits < 1) stbi__grow_buffer_unsafe(j);
   k = j->code_buffer;
   j->code_buffer <<= 1;
   --j->code_bits;
   return k & 0x80000000;
}

// given a value that's at position X in the zigzag stream,
// where does it appear in the 8x8 matrix coded as row-major?
static stbi_uc stbi__jpeg_dezigzag[64+15] =
{
    0,  1,  8, 16,  9,  2,  3, 10,
   17, 24, 32, 25, 18, 11,  4,  5,
   12, 19, 26, 33, 40, 48, 41, 34,
   27, 20, 13,  6,  7, 14, 21, 28,
   35, 42, 49, 56, 57, 50, 43, 36,
   29, 22, 15, 23, 30, 37, 44, 51,
   58, 59, 52, 45, 38, 31, 39, 46,
   53, 60, 61, 54, 47, 55, 62, 63,
   // let corrupt input sample past end
   63, 63, 63, 63, 63, 63, 63, 63,
   63, 63, 63, 63, 63, 63, 63
};

// decode one 64-entry block--
static int stbi__jpeg_decode_block(stbi__jpeg *j, short data[64], stbi__huffman *hdc, stbi__huffman *hac, stbi__int16 *fac, int b, stbi_uc *dequant)
{
   int diff,dc,k;
   int t;

   if (j->code_bits < 16) stbi__grow_buffer_unsafe(j);
   t = stbi__jpeg_huff_decode(j, hdc);
   if (t < 0) return stbi__err("bad huffman code","Corrupt JPEG");

   // 0 all the ac values now so we can do it 32-bits at a time
   memset(data,0,64*sizeof(data[0]));

   diff = t ? stbi__extend_receive(j, t) : 0;
   dc = j->img_comp[b].dc_pred + diff;
   j->img_comp[b].dc_pred = dc;
   data[0] = (short) (dc * dequant[0]);

   // decode AC components, see JPEG spec
   k = 1;
   do {
      unsigned int zig;
      int c,r,s;
      if (j->code_bits < 16) stbi__grow_buffer_unsafe(j);
      c = (j->code_buffer >> (32 - FAST_BITS)) & ((1 << FAST_BITS)-1);
      r = fac[c];
      if (r) { // fast-AC path
         k += (r >> 4) & 15; // run
         s = r & 15; // combined length
         j->code_buffer <<= s;
         j->code_bits -= s;
         // decode into unzigzag'd location
         zig = stbi__jpeg_dezigzag[k++];
         data[zig] = (short) ((r >> 8) * dequant[zig]);
      } else {
         int rs = stbi__jpeg_huff_decode(j, hac);
         if (rs < 0) return stbi__err("bad huffman code","Corrupt JPEG");
         s = rs & 15;
         r = rs >> 4;
         if (s == 0) {
            if (rs != 0xf0) break; // end block
            k += 16;
         } else {
            k += r;
            // decode into unzigzag'd location
            zig = stbi__jpeg_dezigzag[k++];
            data[zig] = (short) (stbi__extend_receive(j,s) * dequant[zig]);
         }
      }
   } while (k < 64);
   return 1;
}

static int stbi__jpeg_decode_block_prog_dc(stbi__jpeg *j, short data[64], stbi__huffman *hdc, int b)
{
   int diff,dc;
   int t;
   if (j->spec_end != 0) return stbi__err("can't merge dc and ac", "Corrupt JPEG");

   if (j->code_bits < 16) stbi__grow_buffer_unsafe(j);

   if (j->succ_high == 0) {
      // first scan for DC coefficient, must be first
      memset(data,0,64*sizeof(data[0])); // 0 all the ac values now
      t = stbi__jpeg_huff_decode(j, hdc);
      diff = t ? stbi__extend_receive(j, t) : 0;

      dc = j->img_comp[b].dc_pred + diff;
      j->img_comp[b].dc_pred = dc;
      data[0] = (short) (dc << j->succ_low);
   } else {
      // refinement scan for DC coefficient
      if (stbi__jpeg_get_bit(j))
         data[0] += (short) (1 << j->succ_low);
   }
   return 1;
}

// @OPTIMIZE: store non-zigzagged during the decode passes,
// and only de-zigzag when dequantizing
static int stbi__jpeg_decode_block_prog_ac(stbi__jpeg *j, short data[64], stbi__huffman *hac, stbi__int16 *fac)
{
   int k;
   if (j->spec_start == 0) return stbi__err("can't merge dc and ac", "Corrupt JPEG");

   if (j->succ_high == 0) {
      int shift = j->succ_low;

      if (j->eob_run) {
         --j->eob_run;
         return 1;
      }

      k = j->spec_start;
      do {
         unsigned int zig;
         int c,r,s;
         if (j->code_bits < 16) stbi__grow_buffer_unsafe(j);
         c = (j->code_buffer >> (32 - FAST_BITS)) & ((1 << FAST_BITS)-1);
         r = fac[c];
         if (r) { // fast-AC path
            k += (r >> 4) & 15; // run
            s = r & 15; // combined length
            j->code_buffer <<= s;
            j->code_bits -= s;
            zig = stbi__jpeg_dezigzag[k++];
            data[zig] = (short) ((r >> 8) << shift);
         } else {
            int rs = stbi__jpeg_huff_decode(j, hac);
            if (rs < 0) return stbi__err("bad huffman code","Corrupt JPEG");
            s = rs & 15;
            r = rs >> 4;
            if (s == 0) {
               if (r < 15) {
                  j->eob_run = (1 << r);
                  if (r)
                     j->eob_run += stbi__jpeg_get_bits(j, r);
                  --j->eob_run;
                  break;
               }
               k += 16;
            } else {
               k += r;
               zig = stbi__jpeg_dezigzag[k++];
               data[zig] = (short) (stbi__extend_receive(j,s) << shift);
            }
         }
      } while (k <= j->spec_end);
   } else {
      // refinement scan for these AC coefficients

      short bit = (short) (1 << j->succ_low);

      if (j->eob_run) {
         --j->eob_run;
         for (k = j->spec_start; k <= j->spec_end; ++k) {
            short *p = &data[stbi__jpeg_dezigzag[k]];
            if (*p != 0)
               if (stbi__jpeg_get_bit(j))
                  if ((*p & bit)==0) {
                     if (*p > 0)
                        *p += bit;
                     else
                        *p -= bit;
                  }
         }
      } else {
         k = j->spec_start;
         do {
            int r,s;
            int rs = stbi__jpeg_huff_decode(j, hac); // @OPTIMIZE see if we can use the fast path here, advance-by-r is so slow, eh
            if (rs < 0) return stbi__err("bad huffman code","Corrupt JPEG");
            s = rs & 15;
            r = rs >> 4;
            if (s == 0) {
               if (r < 15) {
                  j->eob_run = (1 << r) - 1;
                  if (r)
                     j->eob_run += stbi__jpeg_get_bits(j, r);
                  r = 64; // force end of block
               } else {
                  // r=15 s=0 should write 16 0s, so we just do
                  // a run of 15 0s and then write s (which is 0),
                  // so we don't have to do anything special here
               }
            } else {
               if (s != 1) return stbi__err("bad huffman code", "Corrupt JPEG");
               // sign bit
               if (stbi__jpeg_get_bit(j))
                  s = bit;
               else
                  s = -bit;
            }

            // advance by r
            while (k <= j->spec_end) {
               short *p = &data[stbi__jpeg_dezigzag[k++]];
               if (*p != 0) {
                  if (stbi__jpeg_get_bit(j))
                     if ((*p & bit)==0) {
                        if (*p > 0)
                           *p += bit;
                        else
                           *p -= bit;
                     }
               } else {
                  if (r == 0) {
                     *p = (short) s;
                     break;
                  }
                  --r;
               }
            }
         } while (k <= j->spec_end);
      }
   }
   return 1;
}

// take a -128..127 value and stbi__clamp it and convert to 0..255
stbi_inline static stbi_uc stbi__clamp(int x)
{
   // trick to use a single test to catch both cases
   if ((unsigned int) x > 255) {
      if (x < 0) return 0;
      if (x > 255) return 255;
   }
   return (stbi_uc) x;
}

#define stbi__f2f(x)  ((int) (((x) * 4096 + 0.5)))
#define stbi__fsh(x)  ((x) << 12)

// derived from jidctint -- DCT_ISLOW
#define STBI__IDCT_1D(s0,s1,s2,s3,s4,s5,s6,s7) \
   int t0,t1,t2,t3,p1,p2,p3,p4,p5,x0,x1,x2,x3; \
   p2 = s2;                                    \
   p3 = s6;                                    \
   p1 = (p2+p3) * stbi__f2f(0.5411961f);       \
   t2 = p1 + p3*stbi__f2f(-1.847759065f);      \
   t3 = p1 + p2*stbi__f2f( 0.765366865f);      \
   p2 = s0;                                    \
   p3 = s4;                                    \
   t0 = stbi__fsh(p2+p3);                      \
   t1 = stbi__fsh(p2-p3);                      \
   x0 = t0+t3;                                 \
   x3 = t0-t3;                                 \
   x1 = t1+t2;                                 \
   x2 = t1-t2;                                 \
   t0 = s7;                                    \
   t1 = s5;                                    \
   t2 = s3;                                    \
   t3 = s1;                                    \
   p3 = t0+t2;                                 \
   p4 = t1+t3;                                 \
   p1 = t0+t3;                                 \
   p2 = t1+t2;                                 \
   p5 = (p3+p4)*stbi__f2f( 1.175875602f);      \
   t0 = t0*stbi__f2f( 0.298631336f);           \
   t1 = t1*stbi__f2f( 2.053119869f);           \
   t2 = t2*stbi__f2f( 3.072711026f);           \
   t3 = t3*stbi__f2f( 1.501321110f);           \
   p1 = p5 + p1*stbi__f2f(-0.899976223f);      \
   p2 = p5 + p2*stbi__f2f(-2.562915447f);      \
   p3 = p3*stbi__f2f(-1.961570560f);           \
   p4 = p4*stbi__f2f(-0.390180644f);           \
   t3 += p1+p4;                                \
   t2 += p2+p3;                                \
   t1 += p2+p4;                                \
   t0 += p1+p3;

static void stbi__idct_block(stbi_uc *out, int out_stride, short data[64])
{
   int i,val[64],*v=val;
   stbi_uc *o;
   short *d = data;

   // columns
   for (i=0; i < 8; ++i,++d, ++v) {
      // if all zeroes, shortcut -- this avoids dequantizing 0s and IDCTing
      if (d[ 8]==0 && d[16]==0 && d[24]==0 && d[32]==0
           && d[40]==0 && d[48]==0 && d[56]==0) {
         //    no shortcut                 0     seconds
         //    (1|2|3|4|5|6|7)==0          0     seconds
         //    all separate               -0.047 seconds
         //    1 && 2|3 && 4|5 && 6|7:    -0.047 seconds
         int dcterm = d[0] << 2;
         v[0] = v[8] = v[16] = v[24] = v[32] = v[40] = v[48] = v[56] = dcterm;
      } else {
         STBI__IDCT_1D(d[ 0],d[ 8],d[16],d[24],d[32],d[40],d[48],d[56])
         // constants scaled things up by 1<<12; let's bring them back
         // down, but keep 2 extra bits of precision
         x0 += 512; x1 += 512; x2 += 512; x3 += 512;
         v[ 0] = (x0+t3) >> 10;
         v[56] = (x0-t3) >> 10;
         v[ 8] = (x1+t2) >> 10;
         v[48] = (x1-t2) >> 10;
         v[16] = (x2+t1) >> 10;
         v[40] = (x2-t1) >> 10;
         v[24] = (x3+t0) >> 10;
         v[32] = (x3-t0) >> 10;
      }
   }

   for (i=0, v=val, o=out; i < 8; ++i,v+=8,o+=out_stride) {
      // no fast case since the first 1D IDCT spread components out
      STBI__IDCT_1D(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7])
      // constants scaled things up by 1<<12, plus we had 1<<2 from first
      // loop, plus horizontal and vertical each scale by sqrt(8) so together
      // we've got an extra 1<<3, so 1<<17 total we need to remove.
      // so we want to round that, which means adding 0.5 * 1<<17,
      // aka 65536. Also, we'll end up with -128 to 127 that we want
      // to encode as 0..255 by adding 128, so we'll add that before the shift
      x0 += 65536 + (128<<17);
      x1 += 65536 + (128<<17);
      x2 += 65536 + (128<<17);
      x3 += 65536 + (128<<17);
      // tried computing the shifts into temps, or'ing the temps to see
      // if any were out of range, but that was slower
      o[0] = stbi__clamp((x0+t3) >> 17);
      o[7] = stbi__clamp((x0-t3) >> 17);
      o[1] = stbi__clamp((x1+t2) >> 17);
      o[6] = stbi__clamp((x1-t2) >> 17);
      o[2] = stbi__clamp((x2+t1) >> 17);
      o[5] = stbi__clamp((x2-t1) >> 17);
      o[3] = stbi__clamp((x3+t0) >> 17);
      o[4] = stbi__clamp((x3-t0) >> 17);
   }
}

#ifdef STBI_SSE2
// sse2 integer IDCT. not the fastest possible implementation but it
// produces bit-identical results to the generic C version so it's
// fully "transparent".
static void stbi__idct_simd(stbi_uc *out, int out_stride, short data[64])
{
   // This is constructed to match our regular (generic) integer IDCT exactly.
   __m128i row0, row1, row2, row3, row4, row5, row6, row7;
   __m128i tmp;

   // dot product constant: even elems=x, odd elems=y
   #define dct_const(x,y)  _mm_setr_epi16((x),(y),(x),(y),(x),(y),(x),(y))

   // out(0) = c0[even]*x + c0[odd]*y   (c0, x, y 16-bit, out 32-bit)
   // out(1) = c1[even]*x + c1[odd]*y
   #define dct_rot(out0,out1, x,y,c0,c1) \
      __m128i c0##lo = _mm_unpacklo_epi16((x),(y)); \
      __m128i c0##hi = _mm_unpackhi_epi16((x),(y)); \
      __m128i out0##_l = _mm_madd_epi16(c0##lo, c0); \
      __m128i out0##_h = _mm_madd_epi16(c0##hi, c0); \
      __m128i out1##_l = _mm_madd_epi16(c0##lo, c1); \
      __m128i out1##_h = _mm_madd_epi16(c0##hi, c1)

   // out = in << 12  (in 16-bit, out 32-bit)
   #define dct_widen(out, in) \
      __m128i out##_l = _mm_srai_epi32(_mm_unpacklo_epi16(_mm_setzero_si128(), (in)), 4); \
      __m128i out##_h = _mm_srai_epi32(_mm_unpackhi_epi16(_mm_setzero_si128(), (in)), 4)

   // wide add
   #define dct_wadd(out, a, b) \
      __m128i out##_l = _mm_add_epi32(a##_l, b##_l); \
      __m128i out##_h = _mm_add_epi32(a##_h, b##_h)

   // wide sub
   #define dct_wsub(out, a, b) \
      __m128i out##_l = _mm_sub_epi32(a##_l, b##_l); \
      __m128i out##_h = _mm_sub_epi32(a##_h, b##_h)

   // butterfly a/b, add bias, then shift by "s" and pack
   #define dct_bfly32o(out0, out1, a,b,bias,s) \
      { \
         __m128i abiased_l = _mm_add_epi32(a##_l, bias); \
         __m128i abiased_h = _mm_add_epi32(a##_h, bias); \
         dct_wadd(sum, abiased, b); \
         dct_wsub(dif, abiased, b); \
         out0 = _mm_packs_epi32(_mm_srai_epi32(sum_l, s), _mm_srai_epi32(sum_h, s)); \
         out1 = _mm_packs_epi32(_mm_srai_epi32(dif_l, s), _mm_srai_epi32(dif_h, s)); \
      }

   // 8-bit interleave step (for transposes)
   #define dct_interleave8(a, b) \
      tmp = a; \
      a = _mm_unpacklo_epi8(a, b); \
      b = _mm_unpackhi_epi8(tmp, b)

   // 16-bit interleave step (for transposes)
   #define dct_interleave16(a, b) \
      tmp = a; \
      a = _mm_unpacklo_epi16(a, b); \
      b = _mm_unpackhi_epi16(tmp, b)

   #define dct_pass(bias,shift) \
      { \
         /* even part */ \
         dct_rot(t2e,t3e, row2,row6, rot0_0,rot0_1); \
         __m128i sum04 = _mm_add_epi16(row0, row4); \
         __m128i dif04 = _mm_sub_epi16(row0, row4); \
         dct_widen(t0e, sum04); \
         dct_widen(t1e, dif04); \
         dct_wadd(x0, t0e, t3e); \
         dct_wsub(x3, t0e, t3e); \
         dct_wadd(x1, t1e, t2e); \
         dct_wsub(x2, t1e, t2e); \
         /* odd part */ \
         dct_rot(y0o,y2o, row7,row3, rot2_0,rot2_1); \
         dct_rot(y1o,y3o, row5,row1, rot3_0,rot3_1); \
         __m128i sum17 = _mm_add_epi16(row1, row7); \
         __m128i sum35 = _mm_add_epi16(row3, row5); \
         dct_rot(y4o,y5o, sum17,sum35, rot1_0,rot1_1); \
         dct_wadd(x4, y0o, y4o); \
         dct_wadd(x5, y1o, y5o); \
         dct_wadd(x6, y2o, y5o); \
         dct_wadd(x7, y3o, y4o); \
         dct_bfly32o(row0,row7, x0,x7,bias,shift); \
         dct_bfly32o(row1,row6, x1,x6,bias,shift); \
         dct_bfly32o(row2,row5, x2,x5,bias,shift); \
         dct_bfly32o(row3,row4, x3,x4,bias,shift); \
      }

   __m128i rot0_0 = dct_const(stbi__f2f(0.5411961f), stbi__f2f(0.5411961f) + stbi__f2f(-1.847759065f));
   __m128i rot0_1 = dct_const(stbi__f2f(0.5411961f) + stbi__f2f( 0.765366865f), stbi__f2f(0.5411961f));
   __m128i rot1_0 = dct_const(stbi__f2f(1.175875602f) + stbi__f2f(-0.899976223f), stbi__f2f(1.175875602f));
   __m128i rot1_1 = dct_const(stbi__f2f(1.175875602f), stbi__f2f(1.175875602f) + stbi__f2f(-2.562915447f));
   __m128i rot2_0 = dct_const(stbi__f2f(-1.961570560f) + stbi__f2f( 0.298631336f), stbi__f2f(-1.961570560f));
   __m128i rot2_1 = dct_const(stbi__f2f(-1.961570560f), stbi__f2f(-1.961570560f) + stbi__f2f( 3.072711026f));
   __m128i rot3_0 = dct_const(stbi__f2f(-0.390180644f) + stbi__f2f( 2.053119869f), stbi__f2f(-0.390180644f));
   __m128i rot3_1 = dct_const(stbi__f2f(-0.390180644f), stbi__f2f(-0.390180644f) + stbi__f2f( 1.501321110f));

   // rounding biases in column/row passes, see stbi__idct_block for explanation.
   __m128i bias_0 = _mm_set1_epi32(512);
   __m128i bias_1 = _mm_set1_epi32(65536 + (128<<17));

   // load
   row0 = _mm_load_si128((const __m128i *) (data + 0*8));
   row1 = _mm_load_si128((const __m128i *) (data + 1*8));
   row2 = _mm_load_si128((const __m128i *) (data + 2*8));
   row3 = _mm_load_si128((const __m128i *) (data + 3*8));
   row4 = _mm_load_si128((const __m128i *) (data + 4*8));
   row5 = _mm_load_si128((const __m128i *) (data + 5*8));
   row6 = _mm_load_si128((const __m128i *) (data + 6*8));
   row7 = _mm_load_si128((const __m128i *) (data + 7*8));

   // column pass
   dct_pass(bias_0, 10);

   {
      // 16bit 8x8 transpose pass 1
      dct_interleave16(row0, row4);
      dct_interleave16(row1, row5);
      dct_interleave16(row2, row6);
      dct_interleave16(row3, row7);

      // transpose pass 2
      dct_interleave16(row0, row2);
      dct_interleave16(row1, row3);
      dct_interleave16(row4, row6);
      dct_interleave16(row5, row7);

      // transpose pass 3
      dct_interleave16(row0, row1);
      dct_interleave16(row2, row3);
      dct_interleave16(row4, row5);
      dct_interleave16(row6, row7);
   }

   // row pass
   dct_pass(bias_1, 17);

   {
      // pack
      __m128i p0 = _mm_packus_epi16(row0, row1); // a0a1a2a3...a7b0b1b2b3...b7
      __m128i p1 = _mm_packus_epi16(row2, row3);
      __m128i p2 = _mm_packus_epi16(row4, row5);
      __m128i p3 = _mm_packus_epi16(row6, row7);

      // 8bit 8x8 transpose pass 1
      dct_interleave8(p0, p2); // a0e0a1e1...
      dct_interleave8(p1, p3); // c0g0c1g1...

      // transpose pass 2
      dct_interleave8(p0, p1); // a0c0e0g0...
      dct_interleave8(p2, p3); // b0d0f0h0...

      // transpose pass 3
      dct_interleave8(p0, p2); // a0b0c0d0...
      dct_interleave8(p1, p3); // a4b4c4d4...

      // store
      _mm_storel_epi64((__m128i *) out, p0); out += out_stride;
      _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p0, 0x4e)); out += out_stride;
      _mm_storel_epi64((__m128i *) out, p2); out += out_stride;
      _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p2, 0x4e)); out += out_stride;
      _mm_storel_epi64((__m128i *) out, p1); out += out_stride;
      _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p1, 0x4e)); out += out_stride;
      _mm_storel_epi64((__m128i *) out, p3); out += out_stride;
      _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p3, 0x4e));
   }

#undef dct_const
#undef dct_rot
#undef dct_widen
#undef dct_wadd
#undef dct_wsub
#undef dct_bfly32o
#undef dct_interleave8
#undef dct_interleave16
#undef dct_pass
}

#endif // STBI_SSE2

#ifdef STBI_NEON

// NEON integer IDCT. should produce bit-identical
// results to the generic C version.
static void stbi__idct_simd(stbi_uc *out, int out_stride, short data[64])
{
   int16x8_t row0, row1, row2, row3, row4, row5, row6, row7;

   int16x4_t rot0_0 = vdup_n_s16(stbi__f2f(0.5411961f));
   int16x4_t rot0_1 = vdup_n_s16(stbi__f2f(-1.847759065f));
   int16x4_t rot0_2 = vdup_n_s16(stbi__f2f( 0.765366865f));
   int16x4_t rot1_0 = vdup_n_s16(stbi__f2f( 1.175875602f));
   int16x4_t rot1_1 = vdup_n_s16(stbi__f2f(-0.899976223f));
   int16x4_t rot1_2 = vdup_n_s16(stbi__f2f(-2.562915447f));
   int16x4_t rot2_0 = vdup_n_s16(stbi__f2f(-1.961570560f));
   int16x4_t rot2_1 = vdup_n_s16(stbi__f2f(-0.390180644f));
   int16x4_t rot3_0 = vdup_n_s16(stbi__f2f( 0.298631336f));
   int16x4_t rot3_1 = vdup_n_s16(stbi__f2f( 2.053119869f));
   int16x4_t rot3_2 = vdup_n_s16(stbi__f2f( 3.072711026f));
   int16x4_t rot3_3 = vdup_n_s16(stbi__f2f( 1.501321110f));

#define dct_long_mul(out, inq, coeff) \
   int32x4_t out##_l = vmull_s16(vget_low_s16(inq), coeff); \
   int32x4_t out##_h = vmull_s16(vget_high_s16(inq), coeff)

#define dct_long_mac(out, acc, inq, coeff) \
   int32x4_t out##_l = vmlal_s16(acc##_l, vget_low_s16(inq), coeff); \
   int32x4_t out##_h = vmlal_s16(acc##_h, vget_high_s16(inq), coeff)

#define dct_widen(out, inq) \
   int32x4_t out##_l = vshll_n_s16(vget_low_s16(inq), 12); \
   int32x4_t out##_h = vshll_n_s16(vget_high_s16(inq), 12)

// wide add
#define dct_wadd(out, a, b) \
   int32x4_t out##_l = vaddq_s32(a##_l, b##_l); \
   int32x4_t out##_h = vaddq_s32(a##_h, b##_h)

// wide sub
#define dct_wsub(out, a, b) \
   int32x4_t out##_l = vsubq_s32(a##_l, b##_l); \
   int32x4_t out##_h = vsubq_s32(a##_h, b##_h)

// butterfly a/b, then shift using "shiftop" by "s" and pack
#define dct_bfly32o(out0,out1, a,b,shiftop,s) \
   { \
      dct_wadd(sum, a, b); \
      dct_wsub(dif, a, b); \
      out0 = vcombine_s16(shiftop(sum_l, s), shiftop(sum_h, s)); \
      out1 = vcombine_s16(shiftop(dif_l, s), shiftop(dif_h, s)); \
   }

#define dct_pass(shiftop, shift) \
   { \
      /* even part */ \
      int16x8_t sum26 = vaddq_s16(row2, row6); \
      dct_long_mul(p1e, sum26, rot0_0); \
      dct_long_mac(t2e, p1e, row6, rot0_1); \
      dct_long_mac(t3e, p1e, row2, rot0_2); \
      int16x8_t sum04 = vaddq_s16(row0, row4); \
      int16x8_t dif04 = vsubq_s16(row0, row4); \
      dct_widen(t0e, sum04); \
      dct_widen(t1e, dif04); \
      dct_wadd(x0, t0e, t3e); \
      dct_wsub(x3, t0e, t3e); \
      dct_wadd(x1, t1e, t2e); \
      dct_wsub(x2, t1e, t2e); \
      /* odd part */ \
      int16x8_t sum15 = vaddq_s16(row1, row5); \
      int16x8_t sum17 = vaddq_s16(row1, row7); \
      int16x8_t sum35 = vaddq_s16(row3, row5); \
      int16x8_t sum37 = vaddq_s16(row3, row7); \
      int16x8_t sumodd = vaddq_s16(sum17, sum35); \
      dct_long_mul(p5o, sumodd, rot1_0); \
      dct_long_mac(p1o, p5o, sum17, rot1_1); \
      dct_long_mac(p2o, p5o, sum35, rot1_2); \
      dct_long_mul(p3o, sum37, rot2_0); \
      dct_long_mul(p4o, sum15, rot2_1); \
      dct_wadd(sump13o, p1o, p3o); \
      dct_wadd(sump24o, p2o, p4o); \
      dct_wadd(sump23o, p2o, p3o); \
      dct_wadd(sump14o, p1o, p4o); \
      dct_long_mac(x4, sump13o, row7, rot3_0); \
      dct_long_mac(x5, sump24o, row5, rot3_1); \
      dct_long_mac(x6, sump23o, row3, rot3_2); \
      dct_long_mac(x7, sump14o, row1, rot3_3); \
      dct_bfly32o(row0,row7, x0,x7,shiftop,shift); \
      dct_bfly32o(row1,row6, x1,x6,shiftop,shift); \
      dct_bfly32o(row2,row5, x2,x5,shiftop,shift); \
      dct_bfly32o(row3,row4, x3,x4,shiftop,shift); \
   }

   // load
   row0 = vld1q_s16(data + 0*8);
   row1 = vld1q_s16(data + 1*8);
   row2 = vld1q_s16(data + 2*8);
   row3 = vld1q_s16(data + 3*8);
   row4 = vld1q_s16(data + 4*8);
   row5 = vld1q_s16(data + 5*8);
   row6 = vld1q_s16(data + 6*8);
   row7 = vld1q_s16(data + 7*8);

   // add DC bias
   row0 = vaddq_s16(row0, vsetq_lane_s16(1024, vdupq_n_s16(0), 0));

   // column pass
   dct_pass(vrshrn_n_s32, 10);

   // 16bit 8x8 transpose
   {
// these three map to a single VTRN.16, VTRN.32, and VSWP, respectively.
// whether compilers actually get this is another story, sadly.
#define dct_trn16(x, y) { int16x8x2_t t = vtrnq_s16(x, y); x = t.val[0]; y = t.val[1]; }
#define dct_trn32(x, y) { int32x4x2_t t = vtrnq_s32(vreinterpretq_s32_s16(x), vreinterpretq_s32_s16(y)); x = vreinterpretq_s16_s32(t.val[0]); y = vreinterpretq_s16_s32(t.val[1]); }
#define dct_trn64(x, y) { int16x8_t x0 = x; int16x8_t y0 = y; x = vcombine_s16(vget_low_s16(x0), vget_low_s16(y0)); y = vcombine_s16(vget_high_s16(x0), vget_high_s16(y0)); }

      // pass 1
      dct_trn16(row0, row1); // a0b0a2b2a4b4a6b6
      dct_trn16(row2, row3);
      dct_trn16(row4, row5);
      dct_trn16(row6, row7);

      // pass 2
      dct_trn32(row0, row2); // a0b0c0d0a4b4c4d4
      dct_trn32(row1, row3);
      dct_trn32(row4, row6);
      dct_trn32(row5, row7);

      // pass 3
      dct_trn64(row0, row4); // a0b0c0d0e0f0g0h0
      dct_trn64(row1, row5);
      dct_trn64(row2, row6);
      dct_trn64(row3, row7);

#undef dct_trn16
#undef dct_trn32
#undef dct_trn64
   }

   // row pass
   // vrshrn_n_s32 only supports shifts up to 16, we need
   // 17. so do a non-rounding shift of 16 first then follow
   // up with a rounding shift by 1.
   dct_pass(vshrn_n_s32, 16);

   {
      // pack and round
      uint8x8_t p0 = vqrshrun_n_s16(row0, 1);
      uint8x8_t p1 = vqrshrun_n_s16(row1, 1);
      uint8x8_t p2 = vqrshrun_n_s16(row2, 1);
      uint8x8_t p3 = vqrshrun_n_s16(row3, 1);
      uint8x8_t p4 = vqrshrun_n_s16(row4, 1);
      uint8x8_t p5 = vqrshrun_n_s16(row5, 1);
      uint8x8_t p6 = vqrshrun_n_s16(row6, 1);
      uint8x8_t p7 = vqrshrun_n_s16(row7, 1);

      // again, these can translate into one instruction, but often don't.
#define dct_trn8_8(x, y) { uint8x8x2_t t = vtrn_u8(x, y); x = t.val[0]; y = t.val[1]; }
#define dct_trn8_16(x, y) { uint16x4x2_t t = vtrn_u16(vreinterpret_u16_u8(x), vreinterpret_u16_u8(y)); x = vreinterpret_u8_u16(t.val[0]); y = vreinterpret_u8_u16(t.val[1]); }
#define dct_trn8_32(x, y) { uint32x2x2_t t = vtrn_u32(vreinterpret_u32_u8(x), vreinterpret_u32_u8(y)); x = vreinterpret_u8_u32(t.val[0]); y = vreinterpret_u8_u32(t.val[1]); }

      // sadly can't use interleaved stores here since we only write
      // 8 bytes to each scan line!

      // 8x8 8-bit transpose pass 1
      dct_trn8_8(p0, p1);
      dct_trn8_8(p2, p3);
      dct_trn8_8(p4, p5);
      dct_trn8_8(p6, p7);

      // pass 2
      dct_trn8_16(p0, p2);
      dct_trn8_16(p1, p3);
      dct_trn8_16(p4, p6);
      dct_trn8_16(p5, p7);

      // pass 3
      dct_trn8_32(p0, p4);
      dct_trn8_32(p1, p5);
      dct_trn8_32(p2, p6);
      dct_trn8_32(p3, p7);

      // store
      vst1_u8(out, p0); out += out_stride;
      vst1_u8(out, p1); out += out_stride;
      vst1_u8(out, p2); out += out_stride;
      vst1_u8(out, p3); out += out_stride;
      vst1_u8(out, p4); out += out_stride;
      vst1_u8(out, p5); out += out_stride;
      vst1_u8(out, p6); out += out_stride;
      vst1_u8(out, p7);

#undef dct_trn8_8
#undef dct_trn8_16
#undef dct_trn8_32
   }

#undef dct_long_mul
#undef dct_long_mac
#undef dct_widen
#undef dct_wadd
#undef dct_wsub
#undef dct_bfly32o
#undef dct_pass
}

#endif // STBI_NEON

#define STBI__MARKER_none  0xff
// if there's a pending marker from the entropy stream, return that
// otherwise, fetch from the stream and get a marker. if there's no
// marker, return 0xff, which is never a valid marker value
static stbi_uc stbi__get_marker(stbi__jpeg *j)
{
   stbi_uc x;
   if (j->marker != STBI__MARKER_none) { x = j->marker; j->marker = STBI__MARKER_none; return x; }
   x = stbi__get8(j->s);
   if (x != 0xff) return STBI__MARKER_none;
   while (x == 0xff)
      x = stbi__get8(j->s);
   return x;
}

// in each scan, we'll have scan_n components, and the order
// of the components is specified by order[]
#define STBI__RESTART(x)     ((x) >= 0xd0 && (x) <= 0xd7)

// after a restart interval, stbi__jpeg_reset the entropy decoder and
// the dc prediction
static void stbi__jpeg_reset(stbi__jpeg *j)
{
   j->code_bits = 0;
   j->code_buffer = 0;
   j->nomore = 0;
   j->img_comp[0].dc_pred = j->img_comp[1].dc_pred = j->img_comp[2].dc_pred = 0;
   j->marker = STBI__MARKER_none;
   j->todo = j->restart_interval ? j->restart_interval : 0x7fffffff;
   j->eob_run = 0;
   // no more than 1<<31 MCUs if no restart_interal? that's plenty safe,
   // since we don't even allow 1<<30 pixels
}

static int stbi__parse_entropy_coded_data(stbi__jpeg *z)
{
   stbi__jpeg_reset(z);
   if (!z->progressive) {
      if (z->scan_n == 1) {
         int i,j;
         STBI_SIMD_ALIGN(short, data[64]);
         int n = z->order[0];
         // non-interleaved data, we just need to process one block at a time,
         // in trivial scanline order
         // number of blocks to do just depends on how many actual "pixels" this
         // component has, independent of interleaved MCU blocking and such
         int w = (z->img_comp[n].x+7) >> 3;
         int h = (z->img_comp[n].y+7) >> 3;
         for (j=0; j < h; ++j) {
            for (i=0; i < w; ++i) {
               int ha = z->img_comp[n].ha;
               if (!stbi__jpeg_decode_block(z, data, z->huff_dc+z->img_comp[n].hd, z->huff_ac+ha, z->fast_ac[ha], n, z->dequant[z->img_comp[n].tq])) return 0;
               z->idct_block_kernel(z->img_comp[n].data+z->img_comp[n].w2*j*8+i*8, z->img_comp[n].w2, data);
               // every data block is an MCU, so countdown the restart interval
               if (--z->todo <= 0) {
                  if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                  // if it's NOT a restart, then just bail, so we get corrupt data
                  // rather than no data
                  if (!STBI__RESTART(z->marker)) return 1;
                  stbi__jpeg_reset(z);
               }
            }
         }
         return 1;
      } else { // interleaved
         int i,j,k,x,y;
         STBI_SIMD_ALIGN(short, data[64]);
         for (j=0; j < z->img_mcu_y; ++j) {
            for (i=0; i < z->img_mcu_x; ++i) {
               // scan an interleaved mcu... process scan_n components in order
               for (k=0; k < z->scan_n; ++k) {
                  int n = z->order[k];
                  // scan out an mcu's worth of this component; that's just determined
                  // by the basic H and V specified for the component
                  for (y=0; y < z->img_comp[n].v; ++y) {
                     for (x=0; x < z->img_comp[n].h; ++x) {
                        int x2 = (i*z->img_comp[n].h + x)*8;
                        int y2 = (j*z->img_comp[n].v + y)*8;
                        int ha = z->img_comp[n].ha;
                        if (!stbi__jpeg_decode_block(z, data, z->huff_dc+z->img_comp[n].hd, z->huff_ac+ha, z->fast_ac[ha], n, z->dequant[z->img_comp[n].tq])) return 0;
                        z->idct_block_kernel(z->img_comp[n].data+z->img_comp[n].w2*y2+x2, z->img_comp[n].w2, data);
                     }
                  }
               }
               // after all interleaved components, that's an interleaved MCU,
               // so now count down the restart interval
               if (--z->todo <= 0) {
                  if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                  if (!STBI__RESTART(z->marker)) return 1;
                  stbi__jpeg_reset(z);
               }
            }
         }
         return 1;
      }
   } else {
      if (z->scan_n == 1) {
         int i,j;
         int n = z->order[0];
         // non-interleaved data, we just need to process one block at a time,
         // in trivial scanline order
         // number of blocks to do just depends on how many actual "pixels" this
         // component has, independent of interleaved MCU blocking and such
         int w = (z->img_comp[n].x+7) >> 3;
         int h = (z->img_comp[n].y+7) >> 3;
         for (j=0; j < h; ++j) {
            for (i=0; i < w; ++i) {
               short *data = z->img_comp[n].coeff + 64 * (i + j * z->img_comp[n].coeff_w);
               if (z->spec_start == 0) {
                  if (!stbi__jpeg_decode_block_prog_dc(z, data, &z->huff_dc[z->img_comp[n].hd], n))
                     return 0;
               } else {
                  int ha = z->img_comp[n].ha;
                  if (!stbi__jpeg_decode_block_prog_ac(z, data, &z->huff_ac[ha], z->fast_ac[ha]))
                     return 0;
               }
               // every data block is an MCU, so countdown the restart interval
               if (--z->todo <= 0) {
                  if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                  if (!STBI__RESTART(z->marker)) return 1;
                  stbi__jpeg_reset(z);
               }
            }
         }
         return 1;
      } else { // interleaved
         int i,j,k,x,y;
         for (j=0; j < z->img_mcu_y; ++j) {
            for (i=0; i < z->img_mcu_x; ++i) {
               // scan an interleaved mcu... process scan_n components in order
               for (k=0; k < z->scan_n; ++k) {
                  int n = z->order[k];
                  // scan out an mcu's worth of this component; that's just determined
                  // by the basic H and V specified for the component
                  for (y=0; y < z->img_comp[n].v; ++y) {
                     for (x=0; x < z->img_comp[n].h; ++x) {
                        int x2 = (i*z->img_comp[n].h + x);
                        int y2 = (j*z->img_comp[n].v + y);
                        short *data = z->img_comp[n].coeff + 64 * (x2 + y2 * z->img_comp[n].coeff_w);
                        if (!stbi__jpeg_decode_block_prog_dc(z, data, &z->huff_dc[z->img_comp[n].hd], n))
                           return 0;
                     }
                  }
               }
               // after all interleaved components, that's an interleaved MCU,
               // so now count down the restart interval
               if (--z->todo <= 0) {
                  if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                  if (!STBI__RESTART(z->marker)) return 1;
                  stbi__jpeg_reset(z);
               }
            }
         }
         return 1;
      }
   }
}

static void stbi__jpeg_dequantize(short *data, stbi_uc *dequant)
{
   int i;
   for (i=0; i < 64; ++i)
      data[i] *= dequant[i];
}

static void stbi__jpeg_finish(stbi__jpeg *z)
{
   if (z->progressive) {
      // dequantize and idct the data
      int i,j,n;
      for (n=0; n < z->s->img_n; ++n) {
         int w = (z->img_comp[n].x+7) >> 3;
         int h = (z->img_comp[n].y+7) >> 3;
         for (j=0; j < h; ++j) {
            for (i=0; i < w; ++i) {
               short *data = z->img_comp[n].coeff + 64 * (i + j * z->img_comp[n].coeff_w);
               stbi__jpeg_dequantize(data, z->dequant[z->img_comp[n].tq]);
               z->idct_block_kernel(z->img_comp[n].data+z->img_comp[n].w2*j*8+i*8, z->img_comp[n].w2, data);
            }
         }
      }
   }
}

static int stbi__process_marker(stbi__jpeg *z, int m)
{
   int L;
   switch (m) {
      case STBI__MARKER_none: // no marker found
         return stbi__err("expected marker","Corrupt JPEG");

      case 0xDD: // DRI - specify restart interval
         if (stbi__get16be(z->s) != 4) return stbi__err("bad DRI len","Corrupt JPEG");
         z->restart_interval = stbi__get16be(z->s);
         return 1;

      case 0xDB: // DQT - define quantization table
         L = stbi__get16be(z->s)-2;
         while (L > 0) {
            int q = stbi__get8(z->s);
            int p = q >> 4;
            int t = q & 15,i;
            if (p != 0) return stbi__err("bad DQT type","Corrupt JPEG");
            if (t > 3) return stbi__err("bad DQT table","Corrupt JPEG");
            for (i=0; i < 64; ++i)
               z->dequant[t][stbi__jpeg_dezigzag[i]] = stbi__get8(z->s);
            L -= 65;
         }
         return L==0;

      case 0xC4: // DHT - define huffman table
         L = stbi__get16be(z->s)-2;
         while (L > 0) {
            stbi_uc *v;
            int sizes[16],i,n=0;
            int q = stbi__get8(z->s);
            int tc = q >> 4;
            int th = q & 15;
            if (tc > 1 || th > 3) return stbi__err("bad DHT header","Corrupt JPEG");
            for (i=0; i < 16; ++i) {
               sizes[i] = stbi__get8(z->s);
               n += sizes[i];
            }
            L -= 17;
            if (tc == 0) {
               if (!stbi__build_huffman(z->huff_dc+th, sizes)) return 0;
               v = z->huff_dc[th].values;
            } else {
               if (!stbi__build_huffman(z->huff_ac+th, sizes)) return 0;
               v = z->huff_ac[th].values;
            }
            for (i=0; i < n; ++i)
               v[i] = stbi__get8(z->s);
            if (tc != 0)
               stbi__build_fast_ac(z->fast_ac[th], z->huff_ac + th);
            L -= n;
         }
         return L==0;
   }
   // check for comment block or APP blocks
   if ((m >= 0xE0 && m <= 0xEF) || m == 0xFE) {
      stbi__skip(z->s, stbi__get16be(z->s)-2);
      return 1;
   }
   return 0;
}

// after we see SOS
static int stbi__process_scan_header(stbi__jpeg *z)
{
   int i;
   int Ls = stbi__get16be(z->s);
   z->scan_n = stbi__get8(z->s);
   if (z->scan_n < 1 || z->scan_n > 4 || z->scan_n > (int) z->s->img_n) return stbi__err("bad SOS component count","Corrupt JPEG");
   if (Ls != 6+2*z->scan_n) return stbi__err("bad SOS len","Corrupt JPEG");
   for (i=0; i < z->scan_n; ++i) {
      int id = stbi__get8(z->s), which;
      int q = stbi__get8(z->s);
      for (which = 0; which < z->s->img_n; ++which)
         if (z->img_comp[which].id == id)
            break;
      if (which == z->s->img_n) return 0; // no match
      z->img_comp[which].hd = q >> 4;   if (z->img_comp[which].hd > 3) return stbi__err("bad DC huff","Corrupt JPEG");
      z->img_comp[which].ha = q & 15;   if (z->img_comp[which].ha > 3) return stbi__err("bad AC huff","Corrupt JPEG");
      z->order[i] = which;
   }

   {
      int aa;
      z->spec_start = stbi__get8(z->s);
      z->spec_end   = stbi__get8(z->s); // should be 63, but might be 0
      aa = stbi__get8(z->s);
      z->succ_high = (aa >> 4);
      z->succ_low  = (aa & 15);
      if (z->progressive) {
         if (z->spec_start > 63 || z->spec_end > 63  || z->spec_start > z->spec_end || z->succ_high > 13 || z->succ_low > 13)
            return stbi__err("bad SOS", "Corrupt JPEG");
      } else {
         if (z->spec_start != 0) return stbi__err("bad SOS","Corrupt JPEG");
         if (z->succ_high != 0 || z->succ_low != 0) return stbi__err("bad SOS","Corrupt JPEG");
         z->spec_end = 63;
      }
   }

   return 1;
}

static int stbi__process_frame_header(stbi__jpeg *z, int scan)
{
   stbi__context *s = z->s;
   int Lf,p,i,q, h_max=1,v_max=1,c;
   Lf = stbi__get16be(s);         if (Lf < 11) return stbi__err("bad SOF len","Corrupt JPEG"); // JPEG
   p  = stbi__get8(s);            if (p != 8) return stbi__err("only 8-bit","JPEG format not supported: 8-bit only"); // JPEG baseline
   s->img_y = stbi__get16be(s);   if (s->img_y == 0) return stbi__err("no header height", "JPEG format not supported: delayed height"); // Legal, but we don't handle it--but neither does IJG
   s->img_x = stbi__get16be(s);   if (s->img_x == 0) return stbi__err("0 width","Corrupt JPEG"); // JPEG requires
   c = stbi__get8(s);
   if (c != 3 && c != 1) return stbi__err("bad component count","Corrupt JPEG");    // JFIF requires
   s->img_n = c;
   for (i=0; i < c; ++i) {
      z->img_comp[i].data = NULL;
      z->img_comp[i].linebuf = NULL;
   }

   if (Lf != 8+3*s->img_n) return stbi__err("bad SOF len","Corrupt JPEG");

   for (i=0; i < s->img_n; ++i) {
      z->img_comp[i].id = stbi__get8(s);
      if (z->img_comp[i].id != i+1)   // JFIF requires
         if (z->img_comp[i].id != i)  // some version of jpegtran outputs non-JFIF-compliant files!
            return stbi__err("bad component ID","Corrupt JPEG");
      q = stbi__get8(s);
      z->img_comp[i].h = (q >> 4);  if (!z->img_comp[i].h || z->img_comp[i].h > 4) return stbi__err("bad H","Corrupt JPEG");
      z->img_comp[i].v = q & 15;    if (!z->img_comp[i].v || z->img_comp[i].v > 4) return stbi__err("bad V","Corrupt JPEG");
      z->img_comp[i].tq = stbi__get8(s);  if (z->img_comp[i].tq > 3) return stbi__err("bad TQ","Corrupt JPEG");
   }

   if (scan != STBI__SCAN_load) return 1;

   if ((1 << 30) / s->img_x / s->img_n < s->img_y) return stbi__err("too large", "Image too large to decode");

   for (i=0; i < s->img_n; ++i) {
      if (z->img_comp[i].h > h_max) h_max = z->img_comp[i].h;
      if (z->img_comp[i].v > v_max) v_max = z->img_comp[i].v;
   }

   // compute interleaved mcu info
   z->img_h_max = h_max;
   z->img_v_max = v_max;
   z->img_mcu_w = h_max * 8;
   z->img_mcu_h = v_max * 8;
   z->img_mcu_x = (s->img_x + z->img_mcu_w-1) / z->img_mcu_w;
   z->img_mcu_y = (s->img_y + z->img_mcu_h-1) / z->img_mcu_h;

   for (i=0; i < s->img_n; ++i) {
      // number of effective pixels (e.g. for non-interleaved MCU)
      z->img_comp[i].x = (s->img_x * z->img_comp[i].h + h_max-1) / h_max;
      z->img_comp[i].y = (s->img_y * z->img_comp[i].v + v_max-1) / v_max;
      // to simplify generation, we'll allocate enough memory to decode
      // the bogus oversized data from using interleaved MCUs and their
      // big blocks (e.g. a 16x16 iMCU on an image of width 33); we won't
      // discard the extra data until colorspace conversion
      z->img_comp[i].w2 = z->img_mcu_x * z->img_comp[i].h * 8;
      z->img_comp[i].h2 = z->img_mcu_y * z->img_comp[i].v * 8;
      z->img_comp[i].raw_data = stbi__malloc(z->img_comp[i].w2 * z->img_comp[i].h2+15);

      if (z->img_comp[i].raw_data == NULL) {
         for(--i; i >= 0; --i) {
            STBI_FREE(z->img_comp[i].raw_data);
            z->img_comp[i].data = NULL;
         }
         return stbi__err("outofmem", "Out of memory");
      }
      // align blocks for idct using mmx/sse
      z->img_comp[i].data = (stbi_uc*) (((size_t) z->img_comp[i].raw_data + 15) & ~15);
      z->img_comp[i].linebuf = NULL;
      if (z->progressive) {
         z->img_comp[i].coeff_w = (z->img_comp[i].w2 + 7) >> 3;
         z->img_comp[i].coeff_h = (z->img_comp[i].h2 + 7) >> 3;
         z->img_comp[i].raw_coeff = STBI_MALLOC(z->img_comp[i].coeff_w * z->img_comp[i].coeff_h * 64 * sizeof(short) + 15);
         z->img_comp[i].coeff = (short*) (((size_t) z->img_comp[i].raw_coeff + 15) & ~15);
      } else {
         z->img_comp[i].coeff = 0;
         z->img_comp[i].raw_coeff = 0;
      }
   }

   return 1;
}

// use comparisons since in some cases we handle more than one case (e.g. SOF)
#define stbi__DNL(x)         ((x) == 0xdc)
#define stbi__SOI(x)         ((x) == 0xd8)
#define stbi__EOI(x)         ((x) == 0xd9)
#define stbi__SOF(x)         ((x) == 0xc0 || (x) == 0xc1 || (x) == 0xc2)
#define stbi__SOS(x)         ((x) == 0xda)

#define stbi__SOF_progressive(x)   ((x) == 0xc2)

static int stbi__decode_jpeg_header(stbi__jpeg *z, int scan)
{
   int m;
   z->marker = STBI__MARKER_none; // initialize cached marker to empty
   m = stbi__get_marker(z);
   if (!stbi__SOI(m)) return stbi__err("no SOI","Corrupt JPEG");
   if (scan == STBI__SCAN_type) return 1;
   m = stbi__get_marker(z);
   while (!stbi__SOF(m)) {
      if (!stbi__process_marker(z,m)) return 0;
      m = stbi__get_marker(z);
      while (m == STBI__MARKER_none) {
         // some files have extra padding after their blocks, so ok, we'll scan
         if (stbi__at_eof(z->s)) return stbi__err("no SOF", "Corrupt JPEG");
         m = stbi__get_marker(z);
      }
   }
   z->progressive = stbi__SOF_progressive(m);
   if (!stbi__process_frame_header(z, scan)) return 0;
   return 1;
}

// decode image to YCbCr format
static int stbi__decode_jpeg_image(stbi__jpeg *j)
{
   int m;
   for (m = 0; m < 4; m++) {
      j->img_comp[m].raw_data = NULL;
      j->img_comp[m].raw_coeff = NULL;
   }
   j->restart_interval = 0;
   if (!stbi__decode_jpeg_header(j, STBI__SCAN_load)) return 0;
   m = stbi__get_marker(j);
   while (!stbi__EOI(m)) {
      if (stbi__SOS(m)) {
         if (!stbi__process_scan_header(j)) return 0;
         if (!stbi__parse_entropy_coded_data(j)) return 0;
         if (j->marker == STBI__MARKER_none ) {
            // handle 0s at the end of image data from IP Kamera 9060
            while (!stbi__at_eof(j->s)) {
               int x = stbi__get8(j->s);
               if (x == 255) {
                  j->marker = stbi__get8(j->s);
                  break;
               } else if (x != 0) {
                  return stbi__err("junk before marker", "Corrupt JPEG");
               }
            }
            // if we reach eof without hitting a marker, stbi__get_marker() below will fail and we'll eventually return 0
         }
      } else {
         if (!stbi__process_marker(j, m)) return 0;
      }
      m = stbi__get_marker(j);
   }
   if (j->progressive)
      stbi__jpeg_finish(j);
   return 1;
}

// static jfif-centered resampling (across block boundaries)

typedef stbi_uc *(*resample_row_func)(stbi_uc *out, stbi_uc *in0, stbi_uc *in1,
                                    int w, int hs);

#define stbi__div4(x) ((stbi_uc) ((x) >> 2))

static stbi_uc *resample_row_1(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   STBI_NOTUSED(out);
   STBI_NOTUSED(in_far);
   STBI_NOTUSED(w);
   STBI_NOTUSED(hs);
   return in_near;
}

static stbi_uc* stbi__resample_row_v_2(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   // need to generate two samples vertically for every one in input
   int i;
   STBI_NOTUSED(hs);
   for (i=0; i < w; ++i)
      out[i] = stbi__div4(3*in_near[i] + in_far[i] + 2);
   return out;
}

static stbi_uc*  stbi__resample_row_h_2(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   // need to generate two samples horizontally for every one in input
   int i;
   stbi_uc *input = in_near;

   if (w == 1) {
      // if only one sample, can't do any interpolation
      out[0] = out[1] = input[0];
      return out;
   }

   out[0] = input[0];
   out[1] = stbi__div4(input[0]*3 + input[1] + 2);
   for (i=1; i < w-1; ++i) {
      int n = 3*input[i]+2;
      out[i*2+0] = stbi__div4(n+input[i-1]);
      out[i*2+1] = stbi__div4(n+input[i+1]);
   }
   out[i*2+0] = stbi__div4(input[w-2]*3 + input[w-1] + 2);
   out[i*2+1] = input[w-1];

   STBI_NOTUSED(in_far);
   STBI_NOTUSED(hs);

   return out;
}

#define stbi__div16(x) ((stbi_uc) ((x) >> 4))

static stbi_uc *stbi__resample_row_hv_2(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   // need to generate 2x2 samples for every one in input
   int i,t0,t1;
   if (w == 1) {
      out[0] = out[1] = stbi__div4(3*in_near[0] + in_far[0] + 2);
      return out;
   }

   t1 = 3*in_near[0] + in_far[0];
   out[0] = stbi__div4(t1+2);
   for (i=1; i < w; ++i) {
      t0 = t1;
      t1 = 3*in_near[i]+in_far[i];
      out[i*2-1] = stbi__div16(3*t0 + t1 + 8);
      out[i*2  ] = stbi__div16(3*t1 + t0 + 8);
   }
   out[w*2-1] = stbi__div4(t1+2);

   STBI_NOTUSED(hs);

   return out;
}

#if defined(STBI_SSE2) || defined(STBI_NEON)
static stbi_uc *stbi__resample_row_hv_2_simd(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   // need to generate 2x2 samples for every one in input
   int i=0,t0,t1;

   if (w == 1) {
      out[0] = out[1] = stbi__div4(3*in_near[0] + in_far[0] + 2);
      return out;
   }

   t1 = 3*in_near[0] + in_far[0];
   // process groups of 8 pixels for as long as we can.
   // note we can't handle the last pixel in a row in this loop
   // because we need to handle the filter boundary conditions.
   for (; i < ((w-1) & ~7); i += 8) {
#if defined(STBI_SSE2)
      // load and perform the vertical filtering pass
      // this uses 3*x + y = 4*x + (y - x)
      __m128i zero  = _mm_setzero_si128();
      __m128i farb  = _mm_loadl_epi64((__m128i *) (in_far + i));
      __m128i nearb = _mm_loadl_epi64((__m128i *) (in_near + i));
      __m128i farw  = _mm_unpacklo_epi8(farb, zero);
      __m128i nearw = _mm_unpacklo_epi8(nearb, zero);
      __m128i diff  = _mm_sub_epi16(farw, nearw);
      __m128i nears = _mm_slli_epi16(nearw, 2);
      __m128i curr  = _mm_add_epi16(nears, diff); // current row

      // horizontal filter works the same based on shifted vers of current
      // row. "prev" is current row shifted right by 1 pixel; we need to
      // insert the previous pixel value (from t1).
      // "next" is current row shifted left by 1 pixel, with first pixel
      // of next block of 8 pixels added in.
      __m128i prv0 = _mm_slli_si128(curr, 2);
      __m128i nxt0 = _mm_srli_si128(curr, 2);
      __m128i prev = _mm_insert_epi16(prv0, t1, 0);
      __m128i next = _mm_insert_epi16(nxt0, 3*in_near[i+8] + in_far[i+8], 7);

      // horizontal filter, polyphase implementation since it's convenient:
      // even pixels = 3*cur + prev = cur*4 + (prev - cur)
      // odd  pixels = 3*cur + next = cur*4 + (next - cur)
      // note the shared term.
      __m128i bias  = _mm_set1_epi16(8);
      __m128i curs = _mm_slli_epi16(curr, 2);
      __m128i prvd = _mm_sub_epi16(prev, curr);
      __m128i nxtd = _mm_sub_epi16(next, curr);
      __m128i curb = _mm_add_epi16(curs, bias);
      __m128i even = _mm_add_epi16(prvd, curb);
      __m128i odd  = _mm_add_epi16(nxtd, curb);

      // interleave even and odd pixels, then undo scaling.
      __m128i int0 = _mm_unpacklo_epi16(even, odd);
      __m128i int1 = _mm_unpackhi_epi16(even, odd);
      __m128i de0  = _mm_srli_epi16(int0, 4);
      __m128i de1  = _mm_srli_epi16(int1, 4);

      // pack and write output
      __m128i outv = _mm_packus_epi16(de0, de1);
      _mm_storeu_si128((__m128i *) (out + i*2), outv);
#elif defined(STBI_NEON)
      // load and perform the vertical filtering pass
      // this uses 3*x + y = 4*x + (y - x)
      uint8x8_t farb  = vld1_u8(in_far + i);
      uint8x8_t nearb = vld1_u8(in_near + i);
      int16x8_t diff  = vreinterpretq_s16_u16(vsubl_u8(farb, nearb));
      int16x8_t nears = vreinterpretq_s16_u16(vshll_n_u8(nearb, 2));
      int16x8_t curr  = vaddq_s16(nears, diff); // current row

      // horizontal filter works the same based on shifted vers of current
      // row. "prev" is current row shifted right by 1 pixel; we need to
      // insert the previous pixel value (from t1).
      // "next" is current row shifted left by 1 pixel, with first pixel
      // of next block of 8 pixels added in.
      int16x8_t prv0 = vextq_s16(curr, curr, 7);
      int16x8_t nxt0 = vextq_s16(curr, curr, 1);
      int16x8_t prev = vsetq_lane_s16(t1, prv0, 0);
      int16x8_t next = vsetq_lane_s16(3*in_near[i+8] + in_far[i+8], nxt0, 7);

      // horizontal filter, polyphase implementation since it's convenient:
      // even pixels = 3*cur + prev = cur*4 + (prev - cur)
      // odd  pixels = 3*cur + next = cur*4 + (next - cur)
      // note the shared term.
      int16x8_t curs = vshlq_n_s16(curr, 2);
      int16x8_t prvd = vsubq_s16(prev, curr);
      int16x8_t nxtd = vsubq_s16(next, curr);
      int16x8_t even = vaddq_s16(curs, prvd);
      int16x8_t odd  = vaddq_s16(curs, nxtd);

      // undo scaling and round, then store with even/odd phases interleaved
      uint8x8x2_t o;
      o.val[0] = vqrshrun_n_s16(even, 4);
      o.val[1] = vqrshrun_n_s16(odd,  4);
      vst2_u8(out + i*2, o);
#endif

      // "previous" value for next iter
      t1 = 3*in_near[i+7] + in_far[i+7];
   }

   t0 = t1;
   t1 = 3*in_near[i] + in_far[i];
   out[i*2] = stbi__div16(3*t1 + t0 + 8);

   for (++i; i < w; ++i) {
      t0 = t1;
      t1 = 3*in_near[i]+in_far[i];
      out[i*2-1] = stbi__div16(3*t0 + t1 + 8);
      out[i*2  ] = stbi__div16(3*t1 + t0 + 8);
   }
   out[w*2-1] = stbi__div4(t1+2);

   STBI_NOTUSED(hs);

   return out;
}
#endif

static stbi_uc *stbi__resample_row_generic(stbi_uc *out, stbi_uc *in_near, stbi_uc *in_far, int w, int hs)
{
   // resample with nearest-neighbor
   int i,j;
   STBI_NOTUSED(in_far);
   for (i=0; i < w; ++i)
      for (j=0; j < hs; ++j)
         out[i*hs+j] = in_near[i];
   return out;
}

#ifdef STBI_JPEG_OLD
// this is the same YCbCr-to-RGB calculation that stb_image has used
// historically before the algorithm changes in 1.49
#define float2fixed(x)  ((int) ((x) * 65536 + 0.5))
static void stbi__YCbCr_to_RGB_row(stbi_uc *out, const stbi_uc *y, const stbi_uc *pcb, const stbi_uc *pcr, int count, int step)
{
   int i;
   for (i=0; i < count; ++i) {
      int y_fixed = (y[i] << 16) + 32768; // rounding
      int r,g,b;
      int cr = pcr[i] - 128;
      int cb = pcb[i] - 128;
      r = y_fixed + cr*float2fixed(1.40200f);
      g = y_fixed - cr*float2fixed(0.71414f) - cb*float2fixed(0.34414f);
      b = y_fixed                            + cb*float2fixed(1.77200f);
      r >>= 16;
      g >>= 16;
      b >>= 16;
      if ((unsigned) r > 255) { if (r < 0) r = 0; else r = 255; }
      if ((unsigned) g > 255) { if (g < 0) g = 0; else g = 255; }
      if ((unsigned) b > 255) { if (b < 0) b = 0; else b = 255; }
      out[0] = (stbi_uc)r;
      out[1] = (stbi_uc)g;
      out[2] = (stbi_uc)b;
      out[3] = 255;
      out += step;
   }
}
#else
// this is a reduced-precision calculation of YCbCr-to-RGB introduced
// to make sure the code produces the same results in both SIMD and scalar
#define float2fixed(x)  (((int) ((x) * 4096.0f + 0.5f)) << 8)
static void stbi__YCbCr_to_RGB_row(stbi_uc *out, const stbi_uc *y, const stbi_uc *pcb, const stbi_uc *pcr, int count, int step)
{
   int i;
   for (i=0; i < count; ++i) {
      int y_fixed = (y[i] << 20) + (1<<19); // rounding
      int r,g,b;
      int cr = pcr[i] - 128;
      int cb = pcb[i] - 128;
      r = y_fixed +  cr* float2fixed(1.40200f);
      g = y_fixed + (cr*-float2fixed(0.71414f)) + ((cb*-float2fixed(0.34414f)) & 0xffff0000);
      b = y_fixed                               +   cb* float2fixed(1.77200f);
      r >>= 20;
      g >>= 20;
      b >>= 20;
      if ((unsigned) r > 255) { if (r < 0) r = 0; else r = 255; }
      if ((unsigned) g > 255) { if (g < 0) g = 0; else g = 255; }
      if ((unsigned) b > 255) { if (b < 0) b = 0; else b = 255; }
      out[0] = (stbi_uc)r;
      out[1] = (stbi_uc)g;
      out[2] = (stbi_uc)b;
      out[3] = 255;
      out += step;
   }
}
#endif

#if defined(STBI_SSE2) || defined(STBI_NEON)
static void stbi__YCbCr_to_RGB_simd(stbi_uc *out, stbi_uc *y, stbi_uc *pcb, stbi_uc *pcr, int count, int step)
{
   int i = 0;

#ifdef STBI_SSE2
   // step == 3 is pretty ugly on the final interleave, and i'm not convinced
   // it's useful in practice (you wouldn't use it for textures, for example).
   // so just accelerate step == 4 case.
   if (step == 4) {
      // this is a fairly straightforward implementation and not super-optimized.
      __m128i signflip  = _mm_set1_epi8(-0x80);
      __m128i cr_const0 = _mm_set1_epi16(   (short) ( 1.40200f*4096.0f+0.5f));
      __m128i cr_const1 = _mm_set1_epi16( - (short) ( 0.71414f*4096.0f+0.5f));
      __m128i cb_const0 = _mm_set1_epi16( - (short) ( 0.34414f*4096.0f+0.5f));
      __m128i cb_const1 = _mm_set1_epi16(   (short) ( 1.77200f*4096.0f+0.5f));
      __m128i y_bias = _mm_set1_epi8((char) (unsigned char) 128);
      __m128i xw = _mm_set1_epi16(255); // alpha channel

      for (; i+7 < count; i += 8) {
         // load
         __m128i y_bytes = _mm_loadl_epi64((__m128i *) (y+i));
         __m128i cr_bytes = _mm_loadl_epi64((__m128i *) (pcr+i));
         __m128i cb_bytes = _mm_loadl_epi64((__m128i *) (pcb+i));
         __m128i cr_biased = _mm_xor_si128(cr_bytes, signflip); // -128
         __m128i cb_biased = _mm_xor_si128(cb_bytes, signflip); // -128

         // unpack to short (and left-shift cr, cb by 8)
         __m128i yw  = _mm_unpacklo_epi8(y_bias, y_bytes);
         __m128i crw = _mm_unpacklo_epi8(_mm_setzero_si128(), cr_biased);
         __m128i cbw = _mm_unpacklo_epi8(_mm_setzero_si128(), cb_biased);

         // color transform
         __m128i yws = _mm_srli_epi16(yw, 4);
         __m128i cr0 = _mm_mulhi_epi16(cr_const0, crw);
         __m128i cb0 = _mm_mulhi_epi16(cb_const0, cbw);
         __m128i cb1 = _mm_mulhi_epi16(cbw, cb_const1);
         __m128i cr1 = _mm_mulhi_epi16(crw, cr_const1);
         __m128i rws = _mm_add_epi16(cr0, yws);
         __m128i gwt = _mm_add_epi16(cb0, yws);
         __m128i bws = _mm_add_epi16(yws, cb1);
         __m128i gws = _mm_add_epi16(gwt, cr1);

         // descale
         __m128i rw = _mm_srai_epi16(rws, 4);
         __m128i bw = _mm_srai_epi16(bws, 4);
         __m128i gw = _mm_srai_epi16(gws, 4);

         // back to byte, set up for transpose
         __m128i brb = _mm_packus_epi16(rw, bw);
         __m128i gxb = _mm_packus_epi16(gw, xw);

         // transpose to interleave channels
         __m128i t0 = _mm_unpacklo_epi8(brb, gxb);
         __m128i t1 = _mm_unpackhi_epi8(brb, gxb);
         __m128i o0 = _mm_unpacklo_epi16(t0, t1);
         __m128i o1 = _mm_unpackhi_epi16(t0, t1);

         // store
         _mm_storeu_si128((__m128i *) (out + 0), o0);
         _mm_storeu_si128((__m128i *) (out + 16), o1);
         out += 32;
      }
   }
#endif

#ifdef STBI_NEON
   // in this version, step=3 support would be easy to add. but is there demand?
   if (step == 4) {
      // this is a fairly straightforward implementation and not super-optimized.
      uint8x8_t signflip = vdup_n_u8(0x80);
      int16x8_t cr_const0 = vdupq_n_s16(   (short) ( 1.40200f*4096.0f+0.5f));
      int16x8_t cr_const1 = vdupq_n_s16( - (short) ( 0.71414f*4096.0f+0.5f));
      int16x8_t cb_const0 = vdupq_n_s16( - (short) ( 0.34414f*4096.0f+0.5f));
      int16x8_t cb_const1 = vdupq_n_s16(   (short) ( 1.77200f*4096.0f+0.5f));

      for (; i+7 < count; i += 8) {
         // load
         uint8x8_t y_bytes  = vld1_u8(y + i);
         uint8x8_t cr_bytes = vld1_u8(pcr + i);
         uint8x8_t cb_bytes = vld1_u8(pcb + i);
         int8x8_t cr_biased = vreinterpret_s8_u8(vsub_u8(cr_bytes, signflip));
         int8x8_t cb_biased = vreinterpret_s8_u8(vsub_u8(cb_bytes, signflip));

         // expand to s16
         int16x8_t yws = vreinterpretq_s16_u16(vshll_n_u8(y_bytes, 4));
         int16x8_t crw = vshll_n_s8(cr_biased, 7);
         int16x8_t cbw = vshll_n_s8(cb_biased, 7);

         // color transform
         int16x8_t cr0 = vqdmulhq_s16(crw, cr_const0);
         int16x8_t cb0 = vqdmulhq_s16(cbw, cb_const0);
         int16x8_t cr1 = vqdmulhq_s16(crw, cr_const1);
         int16x8_t cb1 = vqdmulhq_s16(cbw, cb_const1);
         int16x8_t rws = vaddq_s16(yws, cr0);
         int16x8_t gws = vaddq_s16(vaddq_s16(yws, cb0), cr1);
         int16x8_t bws = vaddq_s16(yws, cb1);

         // undo scaling, round, convert to byte
         uint8x8x4_t o;
         o.val[0] = vqrshrun_n_s16(rws, 4);
         o.val[1] = vqrshrun_n_s16(gws, 4);
         o.val[2] = vqrshrun_n_s16(bws, 4);
         o.val[3] = vdup_n_u8(255);

         // store, interleaving r/g/b/a
         vst4_u8(out, o);
         out += 8*4;
      }
   }
#endif

   for (; i < count; ++i) {
      int y_fixed = (y[i] << 20) + (1<<19); // rounding
      int r,g,b;
      int cr = pcr[i] - 128;
      int cb = pcb[i] - 128;
      r = y_fixed + cr* float2fixed(1.40200f);
      g = y_fixed + cr*-float2fixed(0.71414f) + ((cb*-float2fixed(0.34414f)) & 0xffff0000);
      b = y_fixed                             +   cb* float2fixed(1.77200f);
      r >>= 20;
      g >>= 20;
      b >>= 20;
      if ((unsigned) r > 255) { if (r < 0) r = 0; else r = 255; }
      if ((unsigned) g > 255) { if (g < 0) g = 0; else g = 255; }
      if ((unsigned) b > 255) { if (b < 0) b = 0; else b = 255; }
      out[0] = (stbi_uc)r;
      out[1] = (stbi_uc)g;
      out[2] = (stbi_uc)b;
      out[3] = 255;
      out += step;
   }
}
#endif

// set up the kernels
static void stbi__setup_jpeg(stbi__jpeg *j)
{
   j->idct_block_kernel = stbi__idct_block;
   j->YCbCr_to_RGB_kernel = stbi__YCbCr_to_RGB_row;
   j->resample_row_hv_2_kernel = stbi__resample_row_hv_2;

#ifdef STBI_SSE2
   if (stbi__sse2_available()) {
      j->idct_block_kernel = stbi__idct_simd;
      #ifndef STBI_JPEG_OLD
      j->YCbCr_to_RGB_kernel = stbi__YCbCr_to_RGB_simd;
      #endif
      j->resample_row_hv_2_kernel = stbi__resample_row_hv_2_simd;
   }
#endif

#ifdef STBI_NEON
   j->idct_block_kernel = stbi__idct_simd;
   #ifndef STBI_JPEG_OLD
   j->YCbCr_to_RGB_kernel = stbi__YCbCr_to_RGB_simd;
   #endif
   j->resample_row_hv_2_kernel = stbi__resample_row_hv_2_simd;
#endif
}

// clean up the temporary component buffers
static void stbi__cleanup_jpeg(stbi__jpeg *j)
{
   int i;
   for (i=0; i < j->s->img_n; ++i) {
      if (j->img_comp[i].raw_data) {
         STBI_FREE(j->img_comp[i].raw_data);
         j->img_comp[i].raw_data = NULL;
         j->img_comp[i].data = NULL;
      }
      if (j->img_comp[i].raw_coeff) {
         STBI_FREE(j->img_comp[i].raw_coeff);
         j->img_comp[i].raw_coeff = 0;
         j->img_comp[i].coeff = 0;
      }
      if (j->img_comp[i].linebuf) {
         STBI_FREE(j->img_comp[i].linebuf);
         j->img_comp[i].linebuf = NULL;
      }
   }
}

typedef struct
{
   resample_row_func resample;
   stbi_uc *line0,*line1;
   int hs,vs;   // expansion factor in each axis
   int w_lores; // horizontal pixels pre-expansion
   int ystep;   // how far through vertical expansion we are
   int ypos;    // which pre-expansion row we're on
} stbi__resample;

static stbi_uc *load_jpeg_image(stbi__jpeg *z, int *out_x, int *out_y, int *comp, int req_comp)
{
   int n, decode_n;
   z->s->img_n = 0; // make stbi__cleanup_jpeg safe

   // validate req_comp
   if (req_comp < 0 || req_comp > 4) return stbi__errpuc("bad req_comp", "Internal error");

   // load a jpeg image from whichever source, but leave in YCbCr format
   if (!stbi__decode_jpeg_image(z)) { stbi__cleanup_jpeg(z); return NULL; }

   // determine actual number of components to generate
   n = req_comp ? req_comp : z->s->img_n;

   if (z->s->img_n == 3 && n < 3)
      decode_n = 1;
   else
      decode_n = z->s->img_n;

   // resample and color-convert
   {
      int k;
      unsigned int i,j;
      stbi_uc *output;
      stbi_uc *coutput[4];

      stbi__resample res_comp[4];

      for (k=0; k < decode_n; ++k) {
         stbi__resample *r = &res_comp[k];

         // allocate line buffer big enough for upsampling off the edges
         // with upsample factor of 4
         z->img_comp[k].linebuf = (stbi_uc *) stbi__malloc(z->s->img_x + 3);
         if (!z->img_comp[k].linebuf) { stbi__cleanup_jpeg(z); return stbi__errpuc("outofmem", "Out of memory"); }

         r->hs      = z->img_h_max / z->img_comp[k].h;
         r->vs      = z->img_v_max / z->img_comp[k].v;
         r->ystep   = r->vs >> 1;
         r->w_lores = (z->s->img_x + r->hs-1) / r->hs;
         r->ypos    = 0;
         r->line0   = r->line1 = z->img_comp[k].data;

         if      (r->hs == 1 && r->vs == 1) r->resample = resample_row_1;
         else if (r->hs == 1 && r->vs == 2) r->resample = stbi__resample_row_v_2;
         else if (r->hs == 2 && r->vs == 1) r->resample = stbi__resample_row_h_2;
         else if (r->hs == 2 && r->vs == 2) r->resample = z->resample_row_hv_2_kernel;
         else                               r->resample = stbi__resample_row_generic;
      }

      // can't error after this so, this is safe
      output = (stbi_uc *) stbi__malloc(n * z->s->img_x * z->s->img_y + 1);
      if (!output) { stbi__cleanup_jpeg(z); return stbi__errpuc("outofmem", "Out of memory"); }

      // now go ahead and resample
      for (j=0; j < z->s->img_y; ++j) {
         stbi_uc *out = output + n * z->s->img_x * j;
         for (k=0; k < decode_n; ++k) {
            stbi__resample *r = &res_comp[k];
            int y_bot = r->ystep >= (r->vs >> 1);
            coutput[k] = r->resample(z->img_comp[k].linebuf,
                                     y_bot ? r->line1 : r->line0,
                                     y_bot ? r->line0 : r->line1,
                                     r->w_lores, r->hs);
            if (++r->ystep >= r->vs) {
               r->ystep = 0;
               r->line0 = r->line1;
               if (++r->ypos < z->img_comp[k].y)
                  r->line1 += z->img_comp[k].w2;
            }
         }
         if (n >= 3) {
            stbi_uc *y = coutput[0];
            if (z->s->img_n == 3) {
               z->YCbCr_to_RGB_kernel(out, y, coutput[1], coutput[2], z->s->img_x, n);
            } else
               for (i=0; i < z->s->img_x; ++i) {
                  out[0] = out[1] = out[2] = y[i];
                  out[3] = 255; // not used if n==3
                  out += n;
               }
         } else {
            stbi_uc *y = coutput[0];
            if (n == 1)
               for (i=0; i < z->s->img_x; ++i) out[i] = y[i];
            else
               for (i=0; i < z->s->img_x; ++i) *out++ = y[i], *out++ = 255;
         }
      }
      stbi__cleanup_jpeg(z);
      *out_x = z->s->img_x;
      *out_y = z->s->img_y;
      if (comp) *comp  = z->s->img_n; // report original components, not output
      return output;
   }
}

static unsigned char *stbi__jpeg_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   stbi__jpeg j;
   j.s = s;
   stbi__setup_jpeg(&j);
   return load_jpeg_image(&j, x,y,comp,req_comp);
}

static int stbi__jpeg_test(stbi__context *s)
{
   int r;
   stbi__jpeg j;
   j.s = s;
   stbi__setup_jpeg(&j);
   r = stbi__decode_jpeg_header(&j, STBI__SCAN_type);
   stbi__rewind(s);
   return r;
}

static int stbi__jpeg_info_raw(stbi__jpeg *j, int *x, int *y, int *comp)
{
   if (!stbi__decode_jpeg_header(j, STBI__SCAN_header)) {
      stbi__rewind( j->s );
      return 0;
   }
   if (x) *x = j->s->img_x;
   if (y) *y = j->s->img_y;
   if (comp) *comp = j->s->img_n;
   return 1;
}

static int stbi__jpeg_info(stbi__context *s, int *x, int *y, int *comp)
{
   stbi__jpeg j;
   j.s = s;
   return stbi__jpeg_info_raw(&j, x, y, comp);
}
#endif

// public domain zlib decode    v0.2  Sean Barrett 2006-11-18
//    simple implementation
//      - all input must be provided in an upfront buffer
//      - all output is written to a single output buffer (can malloc/realloc)
//    performance
//      - fast huffman

#ifndef STBI_NO_ZLIB

// fast-way is faster to check than jpeg huffman, but slow way is slower
#define STBI__ZFAST_BITS  9 // accelerate all cases in default tables
#define STBI__ZFAST_MASK  ((1 << STBI__ZFAST_BITS) - 1)

// zlib-style huffman encoding
// (jpegs packs from left, zlib from right, so can't share code)
typedef struct
{
   stbi__uint16 fast[1 << STBI__ZFAST_BITS];
   stbi__uint16 firstcode[16];
   int maxcode[17];
   stbi__uint16 firstsymbol[16];
   stbi_uc  size[288];
   stbi__uint16 value[288];
} stbi__zhuffman;

stbi_inline static int stbi__bitreverse16(int n)
{
  n = ((n & 0xAAAA) >>  1) | ((n & 0x5555) << 1);
  n = ((n & 0xCCCC) >>  2) | ((n & 0x3333) << 2);
  n = ((n & 0xF0F0) >>  4) | ((n & 0x0F0F) << 4);
  n = ((n & 0xFF00) >>  8) | ((n & 0x00FF) << 8);
  return n;
}

stbi_inline static int stbi__bit_reverse(int v, int bits)
{
   STBI_ASSERT(bits <= 16);
   // to bit reverse n bits, reverse 16 and shift
   // e.g. 11 bits, bit reverse and shift away 5
   return stbi__bitreverse16(v) >> (16-bits);
}

static int stbi__zbuild_huffman(stbi__zhuffman *z, stbi_uc *sizelist, int num)
{
   int i,k=0;
   int code, next_code[16], sizes[17];

   // DEFLATE spec for generating codes
   memset(sizes, 0, sizeof(sizes));
   memset(z->fast, 0, sizeof(z->fast));
   for (i=0; i < num; ++i)
      ++sizes[sizelist[i]];
   sizes[0] = 0;
   for (i=1; i < 16; ++i)
      if (sizes[i] > (1 << i))
         return stbi__err("bad sizes", "Corrupt PNG");
   code = 0;
   for (i=1; i < 16; ++i) {
      next_code[i] = code;
      z->firstcode[i] = (stbi__uint16) code;
      z->firstsymbol[i] = (stbi__uint16) k;
      code = (code + sizes[i]);
      if (sizes[i])
         if (code-1 >= (1 << i)) return stbi__err("bad codelengths","Corrupt PNG");
      z->maxcode[i] = code << (16-i); // preshift for inner loop
      code <<= 1;
      k += sizes[i];
   }
   z->maxcode[16] = 0x10000; // sentinel
   for (i=0; i < num; ++i) {
      int s = sizelist[i];
      if (s) {
         int c = next_code[s] - z->firstcode[s] + z->firstsymbol[s];
         stbi__uint16 fastv = (stbi__uint16) ((s << 9) | i);
         z->size [c] = (stbi_uc     ) s;
         z->value[c] = (stbi__uint16) i;
         if (s <= STBI__ZFAST_BITS) {
            int k = stbi__bit_reverse(next_code[s],s);
            while (k < (1 << STBI__ZFAST_BITS)) {
               z->fast[k] = fastv;
               k += (1 << s);
            }
         }
         ++next_code[s];
      }
   }
   return 1;
}

// zlib-from-memory implementation for PNG reading
//    because PNG allows splitting the zlib stream arbitrarily,
//    and it's annoying structurally to have PNG call ZLIB call PNG,
//    we require PNG read all the IDATs and combine them into a single
//    memory buffer

typedef struct
{
   stbi_uc *zbuffer, *zbuffer_end;
   int num_bits;
   stbi__uint32 code_buffer;

   char *zout;
   char *zout_start;
   char *zout_end;
   int   z_expandable;

   stbi__zhuffman z_length, z_distance;
} stbi__zbuf;

stbi_inline static stbi_uc stbi__zget8(stbi__zbuf *z)
{
   if (z->zbuffer >= z->zbuffer_end) return 0;
   return *z->zbuffer++;
}

static void stbi__fill_bits(stbi__zbuf *z)
{
   do {
      STBI_ASSERT(z->code_buffer < (1U << z->num_bits));
      z->code_buffer |= stbi__zget8(z) << z->num_bits;
      z->num_bits += 8;
   } while (z->num_bits <= 24);
}

stbi_inline static unsigned int stbi__zreceive(stbi__zbuf *z, int n)
{
   unsigned int k;
   if (z->num_bits < n) stbi__fill_bits(z);
   k = z->code_buffer & ((1 << n) - 1);
   z->code_buffer >>= n;
   z->num_bits -= n;
   return k;
}

static int stbi__zhuffman_decode_slowpath(stbi__zbuf *a, stbi__zhuffman *z)
{
   int b,s,k;
   // not resolved by fast table, so compute it the slow way
   // use jpeg approach, which requires MSbits at top
   k = stbi__bit_reverse(a->code_buffer, 16);
   for (s=STBI__ZFAST_BITS+1; ; ++s)
      if (k < z->maxcode[s])
         break;
   if (s == 16) return -1; // invalid code!
   // code size is s, so:
   b = (k >> (16-s)) - z->firstcode[s] + z->firstsymbol[s];
   STBI_ASSERT(z->size[b] == s);
   a->code_buffer >>= s;
   a->num_bits -= s;
   return z->value[b];
}

stbi_inline static int stbi__zhuffman_decode(stbi__zbuf *a, stbi__zhuffman *z)
{
   int b,s;
   if (a->num_bits < 16) stbi__fill_bits(a);
   b = z->fast[a->code_buffer & STBI__ZFAST_MASK];
   if (b) {
      s = b >> 9;
      a->code_buffer >>= s;
      a->num_bits -= s;
      return b & 511;
   }
   return stbi__zhuffman_decode_slowpath(a, z);
}

static int stbi__zexpand(stbi__zbuf *z, char *zout, int n)  // need to make room for n bytes
{
   char *q;
   int cur, limit;
   z->zout = zout;
   if (!z->z_expandable) return stbi__err("output buffer limit","Corrupt PNG");
   cur   = (int) (z->zout     - z->zout_start);
   limit = (int) (z->zout_end - z->zout_start);
   while (cur + n > limit)
      limit *= 2;
   q = (char *) STBI_REALLOC(z->zout_start, limit);
   if (q == NULL) return stbi__err("outofmem", "Out of memory");
   z->zout_start = q;
   z->zout       = q + cur;
   z->zout_end   = q + limit;
   return 1;
}

static int stbi__zlength_base[31] = {
   3,4,5,6,7,8,9,10,11,13,
   15,17,19,23,27,31,35,43,51,59,
   67,83,99,115,131,163,195,227,258,0,0 };

static int stbi__zlength_extra[31]=
{ 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0,0,0 };

static int stbi__zdist_base[32] = { 1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,
257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577,0,0};

static int stbi__zdist_extra[32] =
{ 0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};

static int stbi__parse_huffman_block(stbi__zbuf *a)
{
   char *zout = a->zout;
   for(;;) {
      int z = stbi__zhuffman_decode(a, &a->z_length);
      if (z < 256) {
         if (z < 0) return stbi__err("bad huffman code","Corrupt PNG"); // error in huffman codes
         if (zout >= a->zout_end) {
            if (!stbi__zexpand(a, zout, 1)) return 0;
            zout = a->zout;
         }
         *zout++ = (char) z;
      } else {
         stbi_uc *p;
         int len,dist;
         if (z == 256) {
            a->zout = zout;
            return 1;
         }
         z -= 257;
         len = stbi__zlength_base[z];
         if (stbi__zlength_extra[z]) len += stbi__zreceive(a, stbi__zlength_extra[z]);
         z = stbi__zhuffman_decode(a, &a->z_distance);
         if (z < 0) return stbi__err("bad huffman code","Corrupt PNG");
         dist = stbi__zdist_base[z];
         if (stbi__zdist_extra[z]) dist += stbi__zreceive(a, stbi__zdist_extra[z]);
         if (zout - a->zout_start < dist) return stbi__err("bad dist","Corrupt PNG");
         if (zout + len > a->zout_end) {
            if (!stbi__zexpand(a, zout, len)) return 0;
            zout = a->zout;
         }
         p = (stbi_uc *) (zout - dist);
         if (dist == 1) { // run of one byte; common in images.
            stbi_uc v = *p;
            if (len) { do *zout++ = v; while (--len); }
         } else {
            if (len) { do *zout++ = *p++; while (--len); }
         }
      }
   }
}

static int stbi__compute_huffman_codes(stbi__zbuf *a)
{
   static stbi_uc length_dezigzag[19] = { 16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15 };
   stbi__zhuffman z_codelength;
   stbi_uc lencodes[286+32+137];//padding for maximum single op
   stbi_uc codelength_sizes[19];
   int i,n;

   int hlit  = stbi__zreceive(a,5) + 257;
   int hdist = stbi__zreceive(a,5) + 1;
   int hclen = stbi__zreceive(a,4) + 4;

   memset(codelength_sizes, 0, sizeof(codelength_sizes));
   for (i=0; i < hclen; ++i) {
      int s = stbi__zreceive(a,3);
      codelength_sizes[length_dezigzag[i]] = (stbi_uc) s;
   }
   if (!stbi__zbuild_huffman(&z_codelength, codelength_sizes, 19)) return 0;

   n = 0;
   while (n < hlit + hdist) {
      int c = stbi__zhuffman_decode(a, &z_codelength);
      if (c < 0 || c >= 19) return stbi__err("bad codelengths", "Corrupt PNG");
      if (c < 16)
         lencodes[n++] = (stbi_uc) c;
      else if (c == 16) {
         c = stbi__zreceive(a,2)+3;
         memset(lencodes+n, lencodes[n-1], c);
         n += c;
      } else if (c == 17) {
         c = stbi__zreceive(a,3)+3;
         memset(lencodes+n, 0, c);
         n += c;
      } else {
         STBI_ASSERT(c == 18);
         c = stbi__zreceive(a,7)+11;
         memset(lencodes+n, 0, c);
         n += c;
      }
   }
   if (n != hlit+hdist) return stbi__err("bad codelengths","Corrupt PNG");
   if (!stbi__zbuild_huffman(&a->z_length, lencodes, hlit)) return 0;
   if (!stbi__zbuild_huffman(&a->z_distance, lencodes+hlit, hdist)) return 0;
   return 1;
}

static int stbi__parse_uncomperssed_block(stbi__zbuf *a)
{
   stbi_uc header[4];
   int len,nlen,k;
   if (a->num_bits & 7)
      stbi__zreceive(a, a->num_bits & 7); // discard
   // drain the bit-packed data into header
   k = 0;
   while (a->num_bits > 0) {
      header[k++] = (stbi_uc) (a->code_buffer & 255); // suppress MSVC run-time check
      a->code_buffer >>= 8;
      a->num_bits -= 8;
   }
   STBI_ASSERT(a->num_bits == 0);
   // now fill header the normal way
   while (k < 4)
      header[k++] = stbi__zget8(a);
   len  = header[1] * 256 + header[0];
   nlen = header[3] * 256 + header[2];
   if (nlen != (len ^ 0xffff)) return stbi__err("zlib corrupt","Corrupt PNG");
   if (a->zbuffer + len > a->zbuffer_end) return stbi__err("read past buffer","Corrupt PNG");
   if (a->zout + len > a->zout_end)
      if (!stbi__zexpand(a, a->zout, len)) return 0;
   memcpy(a->zout, a->zbuffer, len);
   a->zbuffer += len;
   a->zout += len;
   return 1;
}

static int stbi__parse_zlib_header(stbi__zbuf *a)
{
   int cmf   = stbi__zget8(a);
   int cm    = cmf & 15;
   /* int cinfo = cmf >> 4; */
   int flg   = stbi__zget8(a);
   if ((cmf*256+flg) % 31 != 0) return stbi__err("bad zlib header","Corrupt PNG"); // zlib spec
   if (flg & 32) return stbi__err("no preset dict","Corrupt PNG"); // preset dictionary not allowed in png
   if (cm != 8) return stbi__err("bad compression","Corrupt PNG"); // DEFLATE required for png
   // window = 1 << (8 + cinfo)... but who cares, we fully buffer output
   return 1;
}

// should statically initialize these for optimal thread safety
static stbi_uc stbi__zdefault_length[288], stbi__zdefault_distance[32];
static void stbi__init_zdefaults(void)
{
   int i;   // use <= to match clearly with spec
   for (i=0; i <= 143; ++i)     stbi__zdefault_length[i]   = 8;
   for (   ; i <= 255; ++i)     stbi__zdefault_length[i]   = 9;
   for (   ; i <= 279; ++i)     stbi__zdefault_length[i]   = 7;
   for (   ; i <= 287; ++i)     stbi__zdefault_length[i]   = 8;

   for (i=0; i <=  31; ++i)     stbi__zdefault_distance[i] = 5;
}

static int stbi__parse_zlib(stbi__zbuf *a, int parse_header)
{
   int final, type;
   if (parse_header)
      if (!stbi__parse_zlib_header(a)) return 0;
   a->num_bits = 0;
   a->code_buffer = 0;
   do {
      final = stbi__zreceive(a,1);
      type = stbi__zreceive(a,2);
      if (type == 0) {
         if (!stbi__parse_uncomperssed_block(a)) return 0;
      } else if (type == 3) {
         return 0;
      } else {
         if (type == 1) {
            // use fixed code lengths
            if (!stbi__zdefault_distance[31]) stbi__init_zdefaults();
            if (!stbi__zbuild_huffman(&a->z_length  , stbi__zdefault_length  , 288)) return 0;
            if (!stbi__zbuild_huffman(&a->z_distance, stbi__zdefault_distance,  32)) return 0;
         } else {
            if (!stbi__compute_huffman_codes(a)) return 0;
         }
         if (!stbi__parse_huffman_block(a)) return 0;
      }
   } while (!final);
   return 1;
}

static int stbi__do_zlib(stbi__zbuf *a, char *obuf, int olen, int exp, int parse_header)
{
   a->zout_start = obuf;
   a->zout       = obuf;
   a->zout_end   = obuf + olen;
   a->z_expandable = exp;

   return stbi__parse_zlib(a, parse_header);
}

STBIDEF char *stbi_zlib_decode_malloc_guesssize(char *buffer, int len, int initial_size, int *outlen)
{
   stbi__zbuf a;
   char *p = (char *) stbi__malloc(initial_size);
   if (p == NULL) return NULL;
   a.zbuffer = (stbi_uc *)buffer;
   a.zbuffer_end = (stbi_uc *)buffer + len;
   if (stbi__do_zlib(&a, p, initial_size, 1, 1)) {
      if (outlen) *outlen = (int) (a.zout - a.zout_start);
      return a.zout_start;
   } else {
      STBI_FREE(a.zout_start);
      return NULL;
   }
}

STBIDEF char *stbi_zlib_decode_malloc(char *buffer, int len, int *outlen)
{
   return stbi_zlib_decode_malloc_guesssize(buffer, len, 16384, outlen);
}

STBIDEF char *stbi_zlib_decode_malloc_guesssize_headerflag(char *buffer, int len, int initial_size, int *outlen, int parse_header)
{
   stbi__zbuf a;
   char *p = (char *) stbi__malloc(initial_size);
   if (p == NULL) return NULL;
   a.zbuffer = (stbi_uc *)buffer;
   a.zbuffer_end = (stbi_uc *)buffer + len;
   if (stbi__do_zlib(&a, p, initial_size, 1, parse_header)) {
      if (outlen) *outlen = (int) (a.zout - a.zout_start);
      return a.zout_start;
   } else {
      STBI_FREE(a.zout_start);
      return NULL;
   }
}

STBIDEF int stbi_zlib_decode_buffer(char *obuffer, int olen, char *ibuffer, int ilen)
{
   stbi__zbuf a;
   a.zbuffer = (stbi_uc *)ibuffer;
   a.zbuffer_end = (stbi_uc *)ibuffer + ilen;
   if (stbi__do_zlib(&a, obuffer, olen, 0, 1))
      return (int) (a.zout - a.zout_start);
   else
      return -1;
}

STBIDEF char *stbi_zlib_decode_noheader_malloc(char *buffer, int len, int *outlen)
{
   stbi__zbuf a;
   char *p = (char *) stbi__malloc(16384);
   if (p == NULL) return NULL;
   a.zbuffer = (stbi_uc *)buffer;
   a.zbuffer_end = (stbi_uc *)buffer+len;
   if (stbi__do_zlib(&a, p, 16384, 1, 0)) {
      if (outlen) *outlen = (int) (a.zout - a.zout_start);
      return a.zout_start;
   } else {
      STBI_FREE(a.zout_start);
      return NULL;
   }
}

STBIDEF int stbi_zlib_decode_noheader_buffer(char *obuffer, int olen, char *ibuffer, int ilen)
{
   stbi__zbuf a;
   a.zbuffer = (stbi_uc *)ibuffer;
   a.zbuffer_end = (stbi_uc *)ibuffer + ilen;
   if (stbi__do_zlib(&a, obuffer, olen, 0, 0))
      return (int) (a.zout - a.zout_start);
   else
      return -1;
}
#endif

// public domain "baseline" PNG decoder   v0.10  Sean Barrett 2006-11-18
//    simple implementation
//      - only 8-bit samples
//      - no CRC checking
//      - allocates lots of intermediate memory
//        - avoids problem of streaming data between subsystems
//        - avoids explicit window management
//    performance
//      - uses stb_zlib, a PD zlib implementation with fast huffman decoding

#ifndef STBI_NO_PNG
typedef struct
{
   stbi__uint32 length;
   stbi__uint32 type;
} stbi__pngchunk;

static stbi__pngchunk stbi__get_chunk_header(stbi__context *s)
{
   stbi__pngchunk c;
   c.length = stbi__get32be(s);
   c.type   = stbi__get32be(s);
   return c;
}

static int stbi__check_png_header(stbi__context *s)
{
   static stbi_uc png_sig[8] = { 137,80,78,71,13,10,26,10 };
   int i;
   for (i=0; i < 8; ++i)
      if (stbi__get8(s) != png_sig[i]) return stbi__err("bad png sig","Not a PNG");
   return 1;
}

typedef struct
{
   stbi__context *s;
   stbi_uc *idata, *expanded, *out;
} stbi__png;


enum {
   STBI__F_none=0,
   STBI__F_sub=1,
   STBI__F_up=2,
   STBI__F_avg=3,
   STBI__F_paeth=4,
   // synthetic filters used for first scanline to avoid needing a dummy row of 0s
   STBI__F_avg_first,
   STBI__F_paeth_first
};

static stbi_uc first_row_filter[5] =
{
   STBI__F_none,
   STBI__F_sub,
   STBI__F_none,
   STBI__F_avg_first,
   STBI__F_paeth_first
};

static int stbi__paeth(int a, int b, int c)
{
   int p = a + b - c;
   int pa = abs(p-a);
   int pb = abs(p-b);
   int pc = abs(p-c);
   if (pa <= pb && pa <= pc) return a;
   if (pb <= pc) return b;
   return c;
}

static stbi_uc stbi__depth_scale_table[9] = { 0, 0xff, 0x55, 0, 0x11, 0,0,0, 0x01 };

// create the png data from post-deflated data
static int stbi__create_png_image_raw(stbi__png *a, stbi_uc *raw, stbi__uint32 raw_len, int out_n, stbi__uint32 x, stbi__uint32 y, int depth, int color)
{
   stbi__context *s = a->s;
   stbi__uint32 i,j,stride = x*out_n;
   stbi__uint32 img_len, img_width_bytes;
   int k;
   int img_n = s->img_n; // copy it into a local for later

   STBI_ASSERT(out_n == s->img_n || out_n == s->img_n+1);
   a->out = (stbi_uc *) stbi__malloc(x * y * out_n); // extra bytes to write off the end into
   if (!a->out) return stbi__err("outofmem", "Out of memory");

   img_width_bytes = (((img_n * x * depth) + 7) >> 3);
   img_len = (img_width_bytes + 1) * y;
   if (s->img_x == x && s->img_y == y) {
      if (raw_len != img_len) return stbi__err("not enough pixels","Corrupt PNG");
   } else { // interlaced:
      if (raw_len < img_len) return stbi__err("not enough pixels","Corrupt PNG");
   }

   for (j=0; j < y; ++j) {
      stbi_uc *cur = a->out + stride*j;
      stbi_uc *prior = cur - stride;
      int filter = *raw++;
      int filter_bytes = img_n;
      int width = x;
      if (filter > 4)
         return stbi__err("invalid filter","Corrupt PNG");

      if (depth < 8) {
         STBI_ASSERT(img_width_bytes <= x);
         cur += x*out_n - img_width_bytes; // store output to the rightmost img_len bytes, so we can decode in place
         filter_bytes = 1;
         width = img_width_bytes;
      }

      // if first row, use special filter that doesn't sample previous row
      if (j == 0) filter = first_row_filter[filter];

      // handle first byte explicitly
      for (k=0; k < filter_bytes; ++k) {
         switch (filter) {
            case STBI__F_none       : cur[k] = raw[k]; break;
            case STBI__F_sub        : cur[k] = raw[k]; break;
            case STBI__F_up         : cur[k] = STBI__BYTECAST(raw[k] + prior[k]); break;
            case STBI__F_avg        : cur[k] = STBI__BYTECAST(raw[k] + (prior[k]>>1)); break;
            case STBI__F_paeth      : cur[k] = STBI__BYTECAST(raw[k] + stbi__paeth(0,prior[k],0)); break;
            case STBI__F_avg_first  : cur[k] = raw[k]; break;
            case STBI__F_paeth_first: cur[k] = raw[k]; break;
         }
      }

      if (depth == 8) {
         if (img_n != out_n)
            cur[img_n] = 255; // first pixel
         raw += img_n;
         cur += out_n;
         prior += out_n;
      } else {
         raw += 1;
         cur += 1;
         prior += 1;
      }

      // this is a little gross, so that we don't switch per-pixel or per-component
      if (depth < 8 || img_n == out_n) {
         int nk = (width - 1)*img_n;
         #define CASE(f) \
             case f:     \
                for (k=0; k < nk; ++k)
         switch (filter) {
            // "none" filter turns into a memcpy here; make that explicit.
            case STBI__F_none:         memcpy(cur, raw, nk); break;
            CASE(STBI__F_sub)          cur[k] = STBI__BYTECAST(raw[k] + cur[k-filter_bytes]); break;
            CASE(STBI__F_up)           cur[k] = STBI__BYTECAST(raw[k] + prior[k]); break;
            CASE(STBI__F_avg)          cur[k] = STBI__BYTECAST(raw[k] + ((prior[k] + cur[k-filter_bytes])>>1)); break;
            CASE(STBI__F_paeth)        cur[k] = STBI__BYTECAST(raw[k] + stbi__paeth(cur[k-filter_bytes],prior[k],prior[k-filter_bytes])); break;
            CASE(STBI__F_avg_first)    cur[k] = STBI__BYTECAST(raw[k] + (cur[k-filter_bytes] >> 1)); break;
            CASE(STBI__F_paeth_first)  cur[k] = STBI__BYTECAST(raw[k] + stbi__paeth(cur[k-filter_bytes],0,0)); break;
         }
         #undef CASE
         raw += nk;
      } else {
         STBI_ASSERT(img_n+1 == out_n);
         #define CASE(f) \
             case f:     \
                for (i=x-1; i >= 1; --i, cur[img_n]=255,raw+=img_n,cur+=out_n,prior+=out_n) \
                   for (k=0; k < img_n; ++k)
         switch (filter) {
            CASE(STBI__F_none)         cur[k] = raw[k]; break;
            CASE(STBI__F_sub)          cur[k] = STBI__BYTECAST(raw[k] + cur[k-out_n]); break;
            CASE(STBI__F_up)           cur[k] = STBI__BYTECAST(raw[k] + prior[k]); break;
            CASE(STBI__F_avg)          cur[k] = STBI__BYTECAST(raw[k] + ((prior[k] + cur[k-out_n])>>1)); break;
            CASE(STBI__F_paeth)        cur[k] = STBI__BYTECAST(raw[k] + stbi__paeth(cur[k-out_n],prior[k],prior[k-out_n])); break;
            CASE(STBI__F_avg_first)    cur[k] = STBI__BYTECAST(raw[k] + (cur[k-out_n] >> 1)); break;
            CASE(STBI__F_paeth_first)  cur[k] = STBI__BYTECAST(raw[k] + stbi__paeth(cur[k-out_n],0,0)); break;
         }
         #undef CASE
      }
   }

   // we make a separate pass to expand bits to pixels; for performance,
   // this could run two scanlines behind the above code, so it won't
   // intefere with filtering but will still be in the cache.
   if (depth < 8) {
      for (j=0; j < y; ++j) {
         stbi_uc *cur = a->out + stride*j;
         stbi_uc *in  = a->out + stride*j + x*out_n - img_width_bytes;
         // unpack 1/2/4-bit into a 8-bit buffer. allows us to keep the common 8-bit path optimal at minimal cost for 1/2/4-bit
         // png guarante byte alignment, if width is not multiple of 8/4/2 we'll decode dummy trailing data that will be skipped in the later loop
         stbi_uc scale = (color == 0) ? stbi__depth_scale_table[depth] : 1; // scale grayscale values to 0..255 range

         // note that the final byte might overshoot and write more data than desired.
         // we can allocate enough data that this never writes out of memory, but it
         // could also overwrite the next scanline. can it overwrite non-empty data
         // on the next scanline? yes, consider 1-pixel-wide scanlines with 1-bit-per-pixel.
         // so we need to explicitly clamp the final ones

         if (depth == 4) {
            for (k=x*img_n; k >= 2; k-=2, ++in) {
               *cur++ = scale * ((*in >> 4)       );
               *cur++ = scale * ((*in     ) & 0x0f);
            }
            if (k > 0) *cur++ = scale * ((*in >> 4)       );
         } else if (depth == 2) {
            for (k=x*img_n; k >= 4; k-=4, ++in) {
               *cur++ = scale * ((*in >> 6)       );
               *cur++ = scale * ((*in >> 4) & 0x03);
               *cur++ = scale * ((*in >> 2) & 0x03);
               *cur++ = scale * ((*in     ) & 0x03);
            }
            if (k > 0) *cur++ = scale * ((*in >> 6)       );
            if (k > 1) *cur++ = scale * ((*in >> 4) & 0x03);
            if (k > 2) *cur++ = scale * ((*in >> 2) & 0x03);
         } else if (depth == 1) {
            for (k=x*img_n; k >= 8; k-=8, ++in) {
               *cur++ = scale * ((*in >> 7)       );
               *cur++ = scale * ((*in >> 6) & 0x01);
               *cur++ = scale * ((*in >> 5) & 0x01);
               *cur++ = scale * ((*in >> 4) & 0x01);
               *cur++ = scale * ((*in >> 3) & 0x01);
               *cur++ = scale * ((*in >> 2) & 0x01);
               *cur++ = scale * ((*in >> 1) & 0x01);
               *cur++ = scale * ((*in     ) & 0x01);
            }
            if (k > 0) *cur++ = scale * ((*in >> 7)       );
            if (k > 1) *cur++ = scale * ((*in >> 6) & 0x01);
            if (k > 2) *cur++ = scale * ((*in >> 5) & 0x01);
            if (k > 3) *cur++ = scale * ((*in >> 4) & 0x01);
            if (k > 4) *cur++ = scale * ((*in >> 3) & 0x01);
            if (k > 5) *cur++ = scale * ((*in >> 2) & 0x01);
            if (k > 6) *cur++ = scale * ((*in >> 1) & 0x01);
         }
         if (img_n != out_n) {
            // insert alpha = 255
            stbi_uc *cur = a->out + stride*j;
            int i;
            if (img_n == 1) {
               for (i=x-1; i >= 0; --i) {
                  cur[i*2+1] = 255;
                  cur[i*2+0] = cur[i];
               }
            } else {
               STBI_ASSERT(img_n == 3);
               for (i=x-1; i >= 0; --i) {
                  cur[i*4+3] = 255;
                  cur[i*4+2] = cur[i*3+2];
                  cur[i*4+1] = cur[i*3+1];
                  cur[i*4+0] = cur[i*3+0];
               }
            }
         }
      }
   }

   return 1;
}

static int stbi__create_png_image(stbi__png *a, stbi_uc *image_data, stbi__uint32 image_data_len, int out_n, int depth, int color, int interlaced)
{
   stbi_uc *final;
   int p;
   if (!interlaced)
      return stbi__create_png_image_raw(a, image_data, image_data_len, out_n, a->s->img_x, a->s->img_y, depth, color);

   // de-interlacing
   final = (stbi_uc *) stbi__malloc(a->s->img_x * a->s->img_y * out_n);
   for (p=0; p < 7; ++p) {
      int xorig[] = { 0,4,0,2,0,1,0 };
      int yorig[] = { 0,0,4,0,2,0,1 };
      int xspc[]  = { 8,8,4,4,2,2,1 };
      int yspc[]  = { 8,8,8,4,4,2,2 };
      int i,j,x,y;
      // pass1_x[4] = 0, pass1_x[5] = 1, pass1_x[12] = 1
      x = (a->s->img_x - xorig[p] + xspc[p]-1) / xspc[p];
      y = (a->s->img_y - yorig[p] + yspc[p]-1) / yspc[p];
      if (x && y) {
         stbi__uint32 img_len = ((((a->s->img_n * x * depth) + 7) >> 3) + 1) * y;
         if (!stbi__create_png_image_raw(a, image_data, image_data_len, out_n, x, y, depth, color)) {
            STBI_FREE(final);
            return 0;
         }
         for (j=0; j < y; ++j) {
            for (i=0; i < x; ++i) {
               int out_y = j*yspc[p]+yorig[p];
               int out_x = i*xspc[p]+xorig[p];
               memcpy(final + out_y*a->s->img_x*out_n + out_x*out_n,
                      a->out + (j*x+i)*out_n, out_n);
            }
         }
         STBI_FREE(a->out);
         image_data += img_len;
         image_data_len -= img_len;
      }
   }
   a->out = final;

   return 1;
}

static int stbi__compute_transparency(stbi__png *z, stbi_uc tc[3], int out_n)
{
   stbi__context *s = z->s;
   stbi__uint32 i, pixel_count = s->img_x * s->img_y;
   stbi_uc *p = z->out;

   // compute color-based transparency, assuming we've
   // already got 255 as the alpha value in the output
   STBI_ASSERT(out_n == 2 || out_n == 4);

   if (out_n == 2) {
      for (i=0; i < pixel_count; ++i) {
         p[1] = (p[0] == tc[0] ? 0 : 255);
         p += 2;
      }
   } else {
      for (i=0; i < pixel_count; ++i) {
         if (p[0] == tc[0] && p[1] == tc[1] && p[2] == tc[2])
            p[3] = 0;
         p += 4;
      }
   }
   return 1;
}

static int stbi__expand_png_palette(stbi__png *a, stbi_uc *palette, int len, int pal_img_n)
{
   stbi__uint32 i, pixel_count = a->s->img_x * a->s->img_y;
   stbi_uc *p, *temp_out, *orig = a->out;

   p = (stbi_uc *) stbi__malloc(pixel_count * pal_img_n);
   if (p == NULL) return stbi__err("outofmem", "Out of memory");

   // between here and free(out) below, exitting would leak
   temp_out = p;

   if (pal_img_n == 3) {
      for (i=0; i < pixel_count; ++i) {
         int n = orig[i]*4;
         p[0] = palette[n  ];
         p[1] = palette[n+1];
         p[2] = palette[n+2];
         p += 3;
      }
   } else {
      for (i=0; i < pixel_count; ++i) {
         int n = orig[i]*4;
         p[0] = palette[n  ];
         p[1] = palette[n+1];
         p[2] = palette[n+2];
         p[3] = palette[n+3];
         p += 4;
      }
   }
   STBI_FREE(a->out);
   a->out = temp_out;

   STBI_NOTUSED(len);

   return 1;
}

static int stbi__unpremultiply_on_load = 0;
static int stbi__de_iphone_flag = 0;

STBIDEF void stbi_set_unpremultiply_on_load(int flag_true_if_should_unpremultiply)
{
   stbi__unpremultiply_on_load = flag_true_if_should_unpremultiply;
}

STBIDEF void stbi_convert_iphone_png_to_rgb(int flag_true_if_should_convert)
{
   stbi__de_iphone_flag = flag_true_if_should_convert;
}

static void stbi__de_iphone(stbi__png *z)
{
   stbi__context *s = z->s;
   stbi__uint32 i, pixel_count = s->img_x * s->img_y;
   stbi_uc *p = z->out;

   if (s->img_out_n == 3) {  // convert bgr to rgb
      for (i=0; i < pixel_count; ++i) {
         stbi_uc t = p[0];
         p[0] = p[2];
         p[2] = t;
         p += 3;
      }
   } else {
      STBI_ASSERT(s->img_out_n == 4);
      if (stbi__unpremultiply_on_load) {
         // convert bgr to rgb and unpremultiply
         for (i=0; i < pixel_count; ++i) {
            stbi_uc a = p[3];
            stbi_uc t = p[0];
            if (a) {
               p[0] = p[2] * 255 / a;
               p[1] = p[1] * 255 / a;
               p[2] =  t   * 255 / a;
            } else {
               p[0] = p[2];
               p[2] = t;
            }
            p += 4;
         }
      } else {
         // convert bgr to rgb
         for (i=0; i < pixel_count; ++i) {
            stbi_uc t = p[0];
            p[0] = p[2];
            p[2] = t;
            p += 4;
         }
      }
   }
}

#define STBI__PNG_TYPE(a,b,c,d)  (((a) << 24) + ((b) << 16) + ((c) << 8) + (d))

static int stbi__parse_png_file(stbi__png *z, int scan, int req_comp)
{
   stbi_uc palette[1024], pal_img_n=0;
   stbi_uc has_trans=0, tc[3];
   stbi__uint32 ioff=0, idata_limit=0, i, pal_len=0;
   int first=1,k,interlace=0, color=0, depth=0, is_iphone=0;
   stbi__context *s = z->s;

   z->expanded = NULL;
   z->idata = NULL;
   z->out = NULL;

   if (!stbi__check_png_header(s)) return 0;

   if (scan == STBI__SCAN_type) return 1;

   for (;;) {
      stbi__pngchunk c = stbi__get_chunk_header(s);
      switch (c.type) {
         case STBI__PNG_TYPE('C','g','B','I'):
            is_iphone = 1;
            stbi__skip(s, c.length);
            break;
         case STBI__PNG_TYPE('I','H','D','R'): {
            int comp,filter;
            if (!first) return stbi__err("multiple IHDR","Corrupt PNG");
            first = 0;
            if (c.length != 13) return stbi__err("bad IHDR len","Corrupt PNG");
            s->img_x = stbi__get32be(s); if (s->img_x > (1 << 24)) return stbi__err("too large","Very large image (corrupt?)");
            s->img_y = stbi__get32be(s); if (s->img_y > (1 << 24)) return stbi__err("too large","Very large image (corrupt?)");
            depth = stbi__get8(s);  if (depth != 1 && depth != 2 && depth != 4 && depth != 8)  return stbi__err("1/2/4/8-bit only","PNG not supported: 1/2/4/8-bit only");
            color = stbi__get8(s);  if (color > 6)         return stbi__err("bad ctype","Corrupt PNG");
            if (color == 3) pal_img_n = 3; else if (color & 1) return stbi__err("bad ctype","Corrupt PNG");
            comp  = stbi__get8(s);  if (comp) return stbi__err("bad comp method","Corrupt PNG");
            filter= stbi__get8(s);  if (filter) return stbi__err("bad filter method","Corrupt PNG");
            interlace = stbi__get8(s); if (interlace>1) return stbi__err("bad interlace method","Corrupt PNG");
            if (!s->img_x || !s->img_y) return stbi__err("0-pixel image","Corrupt PNG");
            if (!pal_img_n) {
               s->img_n = (color & 2 ? 3 : 1) + (color & 4 ? 1 : 0);
               if ((1 << 30) / s->img_x / s->img_n < s->img_y) return stbi__err("too large", "Image too large to decode");
               if (scan == STBI__SCAN_header) return 1;
            } else {
               // if paletted, then pal_n is our final components, and
               // img_n is # components to decompress/filter.
               s->img_n = 1;
               if ((1 << 30) / s->img_x / 4 < s->img_y) return stbi__err("too large","Corrupt PNG");
               // if SCAN_header, have to scan to see if we have a tRNS
            }
            break;
         }

         case STBI__PNG_TYPE('P','L','T','E'):  {
            if (first) return stbi__err("first not IHDR", "Corrupt PNG");
            if (c.length > 256*3) return stbi__err("invalid PLTE","Corrupt PNG");
            pal_len = c.length / 3;
            if (pal_len * 3 != c.length) return stbi__err("invalid PLTE","Corrupt PNG");
            for (i=0; i < pal_len; ++i) {
               palette[i*4+0] = stbi__get8(s);
               palette[i*4+1] = stbi__get8(s);
               palette[i*4+2] = stbi__get8(s);
               palette[i*4+3] = 255;
            }
            break;
         }

         case STBI__PNG_TYPE('t','R','N','S'): {
            if (first) return stbi__err("first not IHDR", "Corrupt PNG");
            if (z->idata) return stbi__err("tRNS after IDAT","Corrupt PNG");
            if (pal_img_n) {
               if (scan == STBI__SCAN_header) { s->img_n = 4; return 1; }
               if (pal_len == 0) return stbi__err("tRNS before PLTE","Corrupt PNG");
               if (c.length > pal_len) return stbi__err("bad tRNS len","Corrupt PNG");
               pal_img_n = 4;
               for (i=0; i < c.length; ++i)
                  palette[i*4+3] = stbi__get8(s);
            } else {
               if (!(s->img_n & 1)) return stbi__err("tRNS with alpha","Corrupt PNG");
               if (c.length != (stbi__uint32) s->img_n*2) return stbi__err("bad tRNS len","Corrupt PNG");
               has_trans = 1;
               for (k=0; k < s->img_n; ++k)
                  tc[k] = (stbi_uc) (stbi__get16be(s) & 255) * stbi__depth_scale_table[depth]; // non 8-bit images will be larger
            }
            break;
         }

         case STBI__PNG_TYPE('I','D','A','T'): {
            if (first) return stbi__err("first not IHDR", "Corrupt PNG");
            if (pal_img_n && !pal_len) return stbi__err("no PLTE","Corrupt PNG");
            if (scan == STBI__SCAN_header) { s->img_n = pal_img_n; return 1; }
            if ((int)(ioff + c.length) < (int)ioff) return 0;
            if (ioff + c.length > idata_limit) {
               stbi_uc *p;
               if (idata_limit == 0) idata_limit = c.length > 4096 ? c.length : 4096;
               while (ioff + c.length > idata_limit)
                  idata_limit *= 2;
               p = (stbi_uc *) STBI_REALLOC(z->idata, idata_limit); if (p == NULL) return stbi__err("outofmem", "Out of memory");
               z->idata = p;
            }
            if (!stbi__getn(s, z->idata+ioff,c.length)) return stbi__err("outofdata","Corrupt PNG");
            ioff += c.length;
            break;
         }

         case STBI__PNG_TYPE('I','E','N','D'): {
            stbi__uint32 raw_len, bpl;
            if (first) return stbi__err("first not IHDR", "Corrupt PNG");
            if (scan != STBI__SCAN_load) return 1;
            if (z->idata == NULL) return stbi__err("no IDAT","Corrupt PNG");
            // initial guess for decoded data size to avoid unnecessary reallocs
            bpl = (s->img_x * depth + 7) / 8; // bytes per line, per component
            raw_len = bpl * s->img_y * s->img_n /* pixels */ + s->img_y /* filter mode per row */;
            z->expanded = (stbi_uc *) stbi_zlib_decode_malloc_guesssize_headerflag((char *) z->idata, ioff, raw_len, (int *) &raw_len, !is_iphone);
            if (z->expanded == NULL) return 0; // zlib should set error
            STBI_FREE(z->idata); z->idata = NULL;
            if ((req_comp == s->img_n+1 && req_comp != 3 && !pal_img_n) || has_trans)
               s->img_out_n = s->img_n+1;
            else
               s->img_out_n = s->img_n;
            if (!stbi__create_png_image(z, z->expanded, raw_len, s->img_out_n, depth, color, interlace)) return 0;
            if (has_trans)
               if (!stbi__compute_transparency(z, tc, s->img_out_n)) return 0;
            if (is_iphone && stbi__de_iphone_flag && s->img_out_n > 2)
               stbi__de_iphone(z);
            if (pal_img_n) {
               // pal_img_n == 3 or 4
               s->img_n = pal_img_n; // record the actual colors we had
               s->img_out_n = pal_img_n;
               if (req_comp >= 3) s->img_out_n = req_comp;
               if (!stbi__expand_png_palette(z, palette, pal_len, s->img_out_n))
                  return 0;
            }
            STBI_FREE(z->expanded); z->expanded = NULL;
            return 1;
         }

         default:
            // if critical, fail
            if (first) return stbi__err("first not IHDR", "Corrupt PNG");
            if ((c.type & (1 << 29)) == 0) {
               #ifndef STBI_NO_FAILURE_STRINGS
               // not threadsafe
               static char invalid_chunk[] = "XXXX PNG chunk not known";
               invalid_chunk[0] = STBI__BYTECAST(c.type >> 24);
               invalid_chunk[1] = STBI__BYTECAST(c.type >> 16);
               invalid_chunk[2] = STBI__BYTECAST(c.type >>  8);
               invalid_chunk[3] = STBI__BYTECAST(c.type >>  0);
               #endif
               return stbi__err(invalid_chunk, "PNG not supported: unknown PNG chunk type");
            }
            stbi__skip(s, c.length);
            break;
      }
      // end of PNG chunk, read and skip CRC
      stbi__get32be(s);
   }
}

static unsigned char *stbi__do_png(stbi__png *p, int *x, int *y, int *n, int req_comp)
{
   unsigned char *result=NULL;
   if (req_comp < 0 || req_comp > 4) return stbi__errpuc("bad req_comp", "Internal error");
   if (stbi__parse_png_file(p, STBI__SCAN_load, req_comp)) {
      result = p->out;
      p->out = NULL;
      if (req_comp && req_comp != p->s->img_out_n) {
         result = stbi__convert_format(result, p->s->img_out_n, req_comp, p->s->img_x, p->s->img_y);
         p->s->img_out_n = req_comp;
         if (result == NULL) return result;
      }
      *x = p->s->img_x;
      *y = p->s->img_y;
      if (n) *n = p->s->img_out_n;
   }
   STBI_FREE(p->out);      p->out      = NULL;
   STBI_FREE(p->expanded); p->expanded = NULL;
   STBI_FREE(p->idata);    p->idata    = NULL;

   return result;
}

static unsigned char *stbi__png_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   stbi__png p;
   p.s = s;
   return stbi__do_png(&p, x,y,comp,req_comp);
}

static int stbi__png_test(stbi__context *s)
{
   int r;
   r = stbi__check_png_header(s);
   stbi__rewind(s);
   return r;
}

static int stbi__png_info_raw(stbi__png *p, int *x, int *y, int *comp)
{
   if (!stbi__parse_png_file(p, STBI__SCAN_header, 0)) {
      stbi__rewind( p->s );
      return 0;
   }
   if (x) *x = p->s->img_x;
   if (y) *y = p->s->img_y;
   if (comp) *comp = p->s->img_n;
   return 1;
}

static int stbi__png_info(stbi__context *s, int *x, int *y, int *comp)
{
   stbi__png p;
   p.s = s;
   return stbi__png_info_raw(&p, x, y, comp);
}
#endif

// Microsoft/Windows BMP image

#ifndef STBI_NO_BMP
static int stbi__bmp_test_raw(stbi__context *s)
{
   int r;
   int sz;
   if (stbi__get8(s) != 'B') return 0;
   if (stbi__get8(s) != 'M') return 0;
   stbi__get32le(s); // discard filesize
   stbi__get16le(s); // discard reserved
   stbi__get16le(s); // discard reserved
   stbi__get32le(s); // discard data offset
   sz = stbi__get32le(s);
   r = (sz == 12 || sz == 40 || sz == 56 || sz == 108 || sz == 124);
   return r;
}

static int stbi__bmp_test(stbi__context *s)
{
   int r = stbi__bmp_test_raw(s);
   stbi__rewind(s);
   return r;
}


// returns 0..31 for the highest set bit
static int stbi__high_bit(unsigned int z)
{
   int n=0;
   if (z == 0) return -1;
   if (z >= 0x10000) n += 16, z >>= 16;
   if (z >= 0x00100) n +=  8, z >>=  8;
   if (z >= 0x00010) n +=  4, z >>=  4;
   if (z >= 0x00004) n +=  2, z >>=  2;
   if (z >= 0x00002) n +=  1, z >>=  1;
   return n;
}

static int stbi__bitcount(unsigned int a)
{
   a = (a & 0x55555555) + ((a >>  1) & 0x55555555); // max 2
   a = (a & 0x33333333) + ((a >>  2) & 0x33333333); // max 4
   a = (a + (a >> 4)) & 0x0f0f0f0f; // max 8 per 4, now 8 bits
   a = (a + (a >> 8)); // max 16 per 8 bits
   a = (a + (a >> 16)); // max 32 per 8 bits
   return a & 0xff;
}

static int stbi__shiftsigned(int v, int shift, int bits)
{
   int result;
   int z=0;

   if (shift < 0) v <<= -shift;
   else v >>= shift;
   result = v;

   z = bits;
   while (z < 8) {
      result += v >> z;
      z += bits;
   }
   return result;
}

static stbi_uc *stbi__bmp_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   stbi_uc *out;
   unsigned int mr=0,mg=0,mb=0,ma=0, fake_a=0;
   stbi_uc pal[256][4];
   int psize=0,i,j,compress=0,width;
   int bpp, flip_vertically, pad, target, offset, hsz;
   if (stbi__get8(s) != 'B' || stbi__get8(s) != 'M') return stbi__errpuc("not BMP", "Corrupt BMP");
   stbi__get32le(s); // discard filesize
   stbi__get16le(s); // discard reserved
   stbi__get16le(s); // discard reserved
   offset = stbi__get32le(s);
   hsz = stbi__get32le(s);
   if (hsz != 12 && hsz != 40 && hsz != 56 && hsz != 108 && hsz != 124) return stbi__errpuc("unknown BMP", "BMP type not supported: unknown");
   if (hsz == 12) {
      s->img_x = stbi__get16le(s);
      s->img_y = stbi__get16le(s);
   } else {
      s->img_x = stbi__get32le(s);
      s->img_y = stbi__get32le(s);
   }
   if (stbi__get16le(s) != 1) return stbi__errpuc("bad BMP", "bad BMP");
   bpp = stbi__get16le(s);
   if (bpp == 1) return stbi__errpuc("monochrome", "BMP type not supported: 1-bit");
   flip_vertically = ((int) s->img_y) > 0;
   s->img_y = abs((int) s->img_y);
   if (hsz == 12) {
      if (bpp < 24)
         psize = (offset - 14 - 24) / 3;
   } else {
      compress = stbi__get32le(s);
      if (compress == 1 || compress == 2) return stbi__errpuc("BMP RLE", "BMP type not supported: RLE");
      stbi__get32le(s); // discard sizeof
      stbi__get32le(s); // discard hres
      stbi__get32le(s); // discard vres
      stbi__get32le(s); // discard colorsused
      stbi__get32le(s); // discard max important
      if (hsz == 40 || hsz == 56) {
         if (hsz == 56) {
            stbi__get32le(s);
            stbi__get32le(s);
            stbi__get32le(s);
            stbi__get32le(s);
         }
         if (bpp == 16 || bpp == 32) {
            mr = mg = mb = 0;
            if (compress == 0) {
               if (bpp == 32) {
                  mr = 0xffu << 16;
                  mg = 0xffu <<  8;
                  mb = 0xffu <<  0;
                  ma = 0xffu << 24;
                  fake_a = 1; // should check for cases like alpha value is all 0 and switch it to 255
                  STBI_NOTUSED(fake_a);
               } else {
                  mr = 31u << 10;
                  mg = 31u <<  5;
                  mb = 31u <<  0;
               }
            } else if (compress == 3) {
               mr = stbi__get32le(s);
               mg = stbi__get32le(s);
               mb = stbi__get32le(s);
               // not documented, but generated by photoshop and handled by mspaint
               if (mr == mg && mg == mb) {
                  // ?!?!?
                  return stbi__errpuc("bad BMP", "bad BMP");
               }
            } else
               return stbi__errpuc("bad BMP", "bad BMP");
         }
      } else {
         STBI_ASSERT(hsz == 108 || hsz == 124);
         mr = stbi__get32le(s);
         mg = stbi__get32le(s);
         mb = stbi__get32le(s);
         ma = stbi__get32le(s);
         stbi__get32le(s); // discard color space
         for (i=0; i < 12; ++i)
            stbi__get32le(s); // discard color space parameters
         if (hsz == 124) {
            stbi__get32le(s); // discard rendering intent
            stbi__get32le(s); // discard offset of profile data
            stbi__get32le(s); // discard size of profile data
            stbi__get32le(s); // discard reserved
         }
      }
      if (bpp < 16)
         psize = (offset - 14 - hsz) >> 2;
   }
   s->img_n = ma ? 4 : 3;
   if (req_comp && req_comp >= 3) // we can directly decode 3 or 4
      target = req_comp;
   else
      target = s->img_n; // if they want monochrome, we'll post-convert
   out = (stbi_uc *) stbi__malloc(target * s->img_x * s->img_y);
   if (!out) return stbi__errpuc("outofmem", "Out of memory");
   if (bpp < 16) {
      int z=0;
      if (psize == 0 || psize > 256) { STBI_FREE(out); return stbi__errpuc("invalid", "Corrupt BMP"); }
      for (i=0; i < psize; ++i) {
         pal[i][2] = stbi__get8(s);
         pal[i][1] = stbi__get8(s);
         pal[i][0] = stbi__get8(s);
         if (hsz != 12) stbi__get8(s);
         pal[i][3] = 255;
      }
      stbi__skip(s, offset - 14 - hsz - psize * (hsz == 12 ? 3 : 4));
      if (bpp == 4) width = (s->img_x + 1) >> 1;
      else if (bpp == 8) width = s->img_x;
      else { STBI_FREE(out); return stbi__errpuc("bad bpp", "Corrupt BMP"); }
      pad = (-width)&3;
      for (j=0; j < (int) s->img_y; ++j) {
         for (i=0; i < (int) s->img_x; i += 2) {
            int v=stbi__get8(s),v2=0;
            if (bpp == 4) {
               v2 = v & 15;
               v >>= 4;
            }
            out[z++] = pal[v][0];
            out[z++] = pal[v][1];
            out[z++] = pal[v][2];
            if (target == 4) out[z++] = 255;
            if (i+1 == (int) s->img_x) break;
            v = (bpp == 8) ? stbi__get8(s) : v2;
            out[z++] = pal[v][0];
            out[z++] = pal[v][1];
            out[z++] = pal[v][2];
            if (target == 4) out[z++] = 255;
         }
         stbi__skip(s, pad);
      }
   } else {
      int rshift=0,gshift=0,bshift=0,ashift=0,rcount=0,gcount=0,bcount=0,acount=0;
      int z = 0;
      int easy=0;
      stbi__skip(s, offset - 14 - hsz);
      if (bpp == 24) width = 3 * s->img_x;
      else if (bpp == 16) width = 2*s->img_x;
      else /* bpp = 32 and pad = 0 */ width=0;
      pad = (-width) & 3;
      if (bpp == 24) {
         easy = 1;
      } else if (bpp == 32) {
         if (mb == 0xff && mg == 0xff00 && mr == 0x00ff0000 && ma == 0xff000000)
            easy = 2;
      }
      if (!easy) {
         if (!mr || !mg || !mb) { STBI_FREE(out); return stbi__errpuc("bad masks", "Corrupt BMP"); }
         // right shift amt to put high bit in position #7
         rshift = stbi__high_bit(mr)-7; rcount = stbi__bitcount(mr);
         gshift = stbi__high_bit(mg)-7; gcount = stbi__bitcount(mg);
         bshift = stbi__high_bit(mb)-7; bcount = stbi__bitcount(mb);
         ashift = stbi__high_bit(ma)-7; acount = stbi__bitcount(ma);
      }
      for (j=0; j < (int) s->img_y; ++j) {
         if (easy) {
            for (i=0; i < (int) s->img_x; ++i) {
               unsigned char a;
               out[z+2] = stbi__get8(s);
               out[z+1] = stbi__get8(s);
               out[z+0] = stbi__get8(s);
               z += 3;
               a = (easy == 2 ? stbi__get8(s) : 255);
               if (target == 4) out[z++] = a;
            }
         } else {
            for (i=0; i < (int) s->img_x; ++i) {
               stbi__uint32 v = (bpp == 16 ? (stbi__uint32) stbi__get16le(s) : stbi__get32le(s));
               int a;
               out[z++] = STBI__BYTECAST(stbi__shiftsigned(v & mr, rshift, rcount));
               out[z++] = STBI__BYTECAST(stbi__shiftsigned(v & mg, gshift, gcount));
               out[z++] = STBI__BYTECAST(stbi__shiftsigned(v & mb, bshift, bcount));
               a = (ma ? stbi__shiftsigned(v & ma, ashift, acount) : 255);
               if (target == 4) out[z++] = STBI__BYTECAST(a);
            }
         }
         stbi__skip(s, pad);
      }
   }
   if (flip_vertically) {
      stbi_uc t;
      for (j=0; j < (int) s->img_y>>1; ++j) {
         stbi_uc *p1 = out +      j     *s->img_x*target;
         stbi_uc *p2 = out + (s->img_y-1-j)*s->img_x*target;
         for (i=0; i < (int) s->img_x*target; ++i) {
            t = p1[i], p1[i] = p2[i], p2[i] = t;
         }
      }
   }

   if (req_comp && req_comp != target) {
      out = stbi__convert_format(out, target, req_comp, s->img_x, s->img_y);
      if (out == NULL) return out; // stbi__convert_format frees input on failure
   }

   *x = s->img_x;
   *y = s->img_y;
   if (comp) *comp = s->img_n;
   return out;
}
#endif

// Targa Truevision - TGA
// by Jonathan Dummer
#ifndef STBI_NO_TGA
static int stbi__tga_info(stbi__context *s, int *x, int *y, int *comp)
{
    int tga_w, tga_h, tga_comp;
    int sz;
    stbi__get8(s);                   // discard Offset
    sz = stbi__get8(s);              // color type
    if( sz > 1 ) {
        stbi__rewind(s);
        return 0;      // only RGB or indexed allowed
    }
    sz = stbi__get8(s);              // image type
    // only RGB or grey allowed, +/- RLE
    if ((sz != 1) && (sz != 2) && (sz != 3) && (sz != 9) && (sz != 10) && (sz != 11)) return 0;
    stbi__skip(s,9);
    tga_w = stbi__get16le(s);
    if( tga_w < 1 ) {
        stbi__rewind(s);
        return 0;   // test width
    }
    tga_h = stbi__get16le(s);
    if( tga_h < 1 ) {
        stbi__rewind(s);
        return 0;   // test height
    }
    sz = stbi__get8(s);               // bits per pixel
    // only RGB or RGBA or grey allowed
    if ((sz != 8) && (sz != 16) && (sz != 24) && (sz != 32)) {
        stbi__rewind(s);
        return 0;
    }
    tga_comp = sz;
    if (x) *x = tga_w;
    if (y) *y = tga_h;
    if (comp) *comp = tga_comp / 8;
    return 1;                   // seems to have passed everything
}

static int stbi__tga_test(stbi__context *s)
{
   int res;
   int sz;
   stbi__get8(s);      //   discard Offset
   sz = stbi__get8(s);   //   color type
   if ( sz > 1 ) return 0;   //   only RGB or indexed allowed
   sz = stbi__get8(s);   //   image type
   if ( (sz != 1) && (sz != 2) && (sz != 3) && (sz != 9) && (sz != 10) && (sz != 11) ) return 0;   //   only RGB or grey allowed, +/- RLE
   stbi__get16be(s);      //   discard palette start
   stbi__get16be(s);      //   discard palette length
   stbi__get8(s);         //   discard bits per palette color entry
   stbi__get16be(s);      //   discard x origin
   stbi__get16be(s);      //   discard y origin
   if ( stbi__get16be(s) < 1 ) return 0;      //   test width
   if ( stbi__get16be(s) < 1 ) return 0;      //   test height
   sz = stbi__get8(s);   //   bits per pixel
   if ( (sz != 8) && (sz != 16) && (sz != 24) && (sz != 32) )
      res = 0;
   else
      res = 1;
   stbi__rewind(s);
   return res;
}

static stbi_uc *stbi__tga_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   //   read in the TGA header stuff
   int tga_offset = stbi__get8(s);
   int tga_indexed = stbi__get8(s);
   int tga_image_type = stbi__get8(s);
   int tga_is_RLE = 0;
   int tga_palette_start = stbi__get16le(s);
   int tga_palette_len = stbi__get16le(s);
   int tga_palette_bits = stbi__get8(s);
   int tga_x_origin = stbi__get16le(s);
   int tga_y_origin = stbi__get16le(s);
   int tga_width = stbi__get16le(s);
   int tga_height = stbi__get16le(s);
   int tga_bits_per_pixel = stbi__get8(s);
   int tga_comp = tga_bits_per_pixel / 8;
   int tga_inverted = stbi__get8(s);
   //   image data
   unsigned char *tga_data;
   unsigned char *tga_palette = NULL;
   int i, j;
   unsigned char raw_data[4];
   int RLE_count = 0;
   int RLE_repeating = 0;
   int read_next_pixel = 1;

   //   do a tiny bit of precessing
   if ( tga_image_type >= 8 )
   {
      tga_image_type -= 8;
      tga_is_RLE = 1;
   }
   /* int tga_alpha_bits = tga_inverted & 15; */
   tga_inverted = 1 - ((tga_inverted >> 5) & 1);

   //   error check
   if ( //(tga_indexed) ||
      (tga_width < 1) || (tga_height < 1) ||
      (tga_image_type < 1) || (tga_image_type > 3) ||
      ((tga_bits_per_pixel != 8) && (tga_bits_per_pixel != 16) &&
      (tga_bits_per_pixel != 24) && (tga_bits_per_pixel != 32))
      )
   {
      return NULL; // we don't report this as a bad TGA because we don't even know if it's TGA
   }

   //   If I'm paletted, then I'll use the number of bits from the palette
   if ( tga_indexed )
   {
      tga_comp = tga_palette_bits / 8;
   }

   //   tga info
   *x = tga_width;
   *y = tga_height;
   if (comp) *comp = tga_comp;

   tga_data = (unsigned char*)stbi__malloc( (size_t)tga_width * tga_height * tga_comp );
   if (!tga_data) return stbi__errpuc("outofmem", "Out of memory");

   // skip to the data's starting position (offset usually = 0)
   stbi__skip(s, tga_offset );

   if ( !tga_indexed && !tga_is_RLE) {
      for (i=0; i < tga_height; ++i) {
         int y = tga_inverted ? tga_height -i - 1 : i;
         stbi_uc *tga_row = tga_data + y*tga_width*tga_comp;
         stbi__getn(s, tga_row, tga_width * tga_comp);
      }
   } else  {
      //   do I need to load a palette?
      if ( tga_indexed)
      {
         //   any data to skip? (offset usually = 0)
         stbi__skip(s, tga_palette_start );
         //   load the palette
         tga_palette = (unsigned char*)stbi__malloc( tga_palette_len * tga_palette_bits / 8 );
         if (!tga_palette) {
            STBI_FREE(tga_data);
            return stbi__errpuc("outofmem", "Out of memory");
         }
         if (!stbi__getn(s, tga_palette, tga_palette_len * tga_palette_bits / 8 )) {
            STBI_FREE(tga_data);
            STBI_FREE(tga_palette);
            return stbi__errpuc("bad palette", "Corrupt TGA");
         }
      }
      //   load the data
      for (i=0; i < tga_width * tga_height; ++i)
      {
         //   if I'm in RLE mode, do I need to get a RLE stbi__pngchunk?
         if ( tga_is_RLE )
         {
            if ( RLE_count == 0 )
            {
               //   yep, get the next byte as a RLE command
               int RLE_cmd = stbi__get8(s);
               RLE_count = 1 + (RLE_cmd & 127);
               RLE_repeating = RLE_cmd >> 7;
               read_next_pixel = 1;
            } else if ( !RLE_repeating )
            {
               read_next_pixel = 1;
            }
         } else
         {
            read_next_pixel = 1;
         }
         //   OK, if I need to read a pixel, do it now
         if ( read_next_pixel )
         {
            //   load however much data we did have
            if ( tga_indexed )
            {
               //   read in 1 byte, then perform the lookup
               int pal_idx = stbi__get8(s);
               if ( pal_idx >= tga_palette_len )
               {
                  //   invalid index
                  pal_idx = 0;
               }
               pal_idx *= tga_bits_per_pixel / 8;
               for (j = 0; j*8 < tga_bits_per_pixel; ++j)
               {
                  raw_data[j] = tga_palette[pal_idx+j];
               }
            } else
            {
               //   read in the data raw
               for (j = 0; j*8 < tga_bits_per_pixel; ++j)
               {
                  raw_data[j] = stbi__get8(s);
               }
            }
            //   clear the reading flag for the next pixel
            read_next_pixel = 0;
         } // end of reading a pixel

         // copy data
         for (j = 0; j < tga_comp; ++j)
           tga_data[i*tga_comp+j] = raw_data[j];

         //   in case we're in RLE mode, keep counting down
         --RLE_count;
      }
      //   do I need to invert the image?
      if ( tga_inverted )
      {
         for (j = 0; j*2 < tga_height; ++j)
         {
            int index1 = j * tga_width * tga_comp;
            int index2 = (tga_height - 1 - j) * tga_width * tga_comp;
            for (i = tga_width * tga_comp; i > 0; --i)
            {
               unsigned char temp = tga_data[index1];
               tga_data[index1] = tga_data[index2];
               tga_data[index2] = temp;
               ++index1;
               ++index2;
            }
         }
      }
      //   clear my palette, if I had one
      if ( tga_palette != NULL )
      {
         STBI_FREE( tga_palette );
      }
   }

   // swap RGB
   if (tga_comp >= 3)
   {
      unsigned char* tga_pixel = tga_data;
      for (i=0; i < tga_width * tga_height; ++i)
      {
         unsigned char temp = tga_pixel[0];
         tga_pixel[0] = tga_pixel[2];
         tga_pixel[2] = temp;
         tga_pixel += tga_comp;
      }
   }

   // convert to target component count
   if (req_comp && req_comp != tga_comp)
      tga_data = stbi__convert_format(tga_data, tga_comp, req_comp, tga_width, tga_height);

   //   the things I do to get rid of an error message, and yet keep
   //   Microsoft's C compilers happy... [8^(
   tga_palette_start = tga_palette_len = tga_palette_bits =
         tga_x_origin = tga_y_origin = 0;
   //   OK, done
   return tga_data;
}
#endif

// *************************************************************************************************
// Photoshop PSD loader -- PD by Thatcher Ulrich, integration by Nicolas Schulz, tweaked by STB

#ifndef STBI_NO_PSD
static int stbi__psd_test(stbi__context *s)
{
   int r = (stbi__get32be(s) == 0x38425053);
   stbi__rewind(s);
   return r;
}

static stbi_uc *stbi__psd_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   int   pixelCount;
   int channelCount, compression;
   int channel, i, count, len;
   int w,h;
   stbi_uc *out;

   // Check identifier
   if (stbi__get32be(s) != 0x38425053)   // "8BPS"
      return stbi__errpuc("not PSD", "Corrupt PSD image");

   // Check file type version.
   if (stbi__get16be(s) != 1)
      return stbi__errpuc("wrong version", "Unsupported version of PSD image");

   // Skip 6 reserved bytes.
   stbi__skip(s, 6 );

   // Read the number of channels (R, G, B, A, etc).
   channelCount = stbi__get16be(s);
   if (channelCount < 0 || channelCount > 16)
      return stbi__errpuc("wrong channel count", "Unsupported number of channels in PSD image");

   // Read the rows and columns of the image.
   h = stbi__get32be(s);
   w = stbi__get32be(s);

   // Make sure the depth is 8 bits.
   if (stbi__get16be(s) != 8)
      return stbi__errpuc("unsupported bit depth", "PSD bit depth is not 8 bit");

   // Make sure the color mode is RGB.
   // Valid options are:
   //   0: Bitmap
   //   1: Grayscale
   //   2: Indexed color
   //   3: RGB color
   //   4: CMYK color
   //   7: Multichannel
   //   8: Duotone
   //   9: Lab color
   if (stbi__get16be(s) != 3)
      return stbi__errpuc("wrong color format", "PSD is not in RGB color format");

   // Skip the Mode Data.  (It's the palette for indexed color; other info for other modes.)
   stbi__skip(s,stbi__get32be(s) );

   // Skip the image resources.  (resolution, pen tool paths, etc)
   stbi__skip(s, stbi__get32be(s) );

   // Skip the reserved data.
   stbi__skip(s, stbi__get32be(s) );

   // Find out if the data is compressed.
   // Known values:
   //   0: no compression
   //   1: RLE compressed
   compression = stbi__get16be(s);
   if (compression > 1)
      return stbi__errpuc("bad compression", "PSD has an unknown compression format");

   // Create the destination image.
   out = (stbi_uc *) stbi__malloc(4 * w*h);
   if (!out) return stbi__errpuc("outofmem", "Out of memory");
   pixelCount = w*h;

   // Initialize the data to zero.
   //memset( out, 0, pixelCount * 4 );

   // Finally, the image data.
   if (compression) {
      // RLE as used by .PSD and .TIFF
      // Loop until you get the number of unpacked bytes you are expecting:
      //     Read the next source byte into n.
      //     If n is between 0 and 127 inclusive, copy the next n+1 bytes literally.
      //     Else if n is between -127 and -1 inclusive, copy the next byte -n+1 times.
      //     Else if n is 128, noop.
      // Endloop

      // The RLE-compressed data is preceeded by a 2-byte data count for each row in the data,
      // which we're going to just skip.
      stbi__skip(s, h * channelCount * 2 );

      // Read the RLE data by channel.
      for (channel = 0; channel < 4; channel++) {
         stbi_uc *p;

         p = out+channel;
         if (channel >= channelCount) {
            // Fill this channel with default data.
            for (i = 0; i < pixelCount; i++, p += 4)
               *p = (channel == 3 ? 255 : 0);
         } else {
            // Read the RLE data.
            count = 0;
            while (count < pixelCount) {
               len = stbi__get8(s);
               if (len == 128) {
                  // No-op.
               } else if (len < 128) {
                  // Copy next len+1 bytes literally.
                  len++;
                  count += len;
                  while (len) {
                     *p = stbi__get8(s);
                     p += 4;
                     len--;
                  }
               } else if (len > 128) {
                  stbi_uc   val;
                  // Next -len+1 bytes in the dest are replicated from next source byte.
                  // (Interpret len as a negative 8-bit int.)
                  len ^= 0x0FF;
                  len += 2;
                  val = stbi__get8(s);
                  count += len;
                  while (len) {
                     *p = val;
                     p += 4;
                     len--;
                  }
               }
            }
         }
      }

   } else {
      // We're at the raw image data.  It's each channel in order (Red, Green, Blue, Alpha, ...)
      // where each channel consists of an 8-bit value for each pixel in the image.

      // Read the data by channel.
      for (channel = 0; channel < 4; channel++) {
         stbi_uc *p;

         p = out + channel;
         if (channel > channelCount) {
            // Fill this channel with default data.
            for (i = 0; i < pixelCount; i++, p += 4)
               *p = channel == 3 ? 255 : 0;
         } else {
            // Read the data.
            for (i = 0; i < pixelCount; i++, p += 4)
               *p = stbi__get8(s);
         }
      }
   }

   if (req_comp && req_comp != 4) {
      out = stbi__convert_format(out, 4, req_comp, w, h);
      if (out == NULL) return out; // stbi__convert_format frees input on failure
   }

   if (comp) *comp = 4;
   *y = h;
   *x = w;

   return out;
}
#endif

// *************************************************************************************************
// Softimage PIC loader
// by Tom Seddon
//
// See http://softimage.wiki.softimage.com/index.php/INFO:_PIC_file_format
// See http://ozviz.wasp.uwa.edu.au/~pbourke/dataformats/softimagepic/

#ifndef STBI_NO_PIC
static int stbi__pic_is4(stbi__context *s,const char *str)
{
   int i;
   for (i=0; i<4; ++i)
      if (stbi__get8(s) != (stbi_uc)str[i])
         return 0;

   return 1;
}

static int stbi__pic_test_core(stbi__context *s)
{
   int i;

   if (!stbi__pic_is4(s,"\x53\x80\xF6\x34"))
      return 0;

   for(i=0;i<84;++i)
      stbi__get8(s);

   if (!stbi__pic_is4(s,"PICT"))
      return 0;

   return 1;
}

typedef struct
{
   stbi_uc size,type,channel;
} stbi__pic_packet;

static stbi_uc *stbi__readval(stbi__context *s, int channel, stbi_uc *dest)
{
   int mask=0x80, i;

   for (i=0; i<4; ++i, mask>>=1) {
      if (channel & mask) {
         if (stbi__at_eof(s)) return stbi__errpuc("bad file","PIC file too short");
         dest[i]=stbi__get8(s);
      }
   }

   return dest;
}

static void stbi__copyval(int channel,stbi_uc *dest,const stbi_uc *src)
{
   int mask=0x80,i;

   for (i=0;i<4; ++i, mask>>=1)
      if (channel&mask)
         dest[i]=src[i];
}

static stbi_uc *stbi__pic_load_core(stbi__context *s,int width,int height,int *comp, stbi_uc *result)
{
   int act_comp=0,num_packets=0,y,chained;
   stbi__pic_packet packets[10];

   // this will (should...) cater for even some bizarre stuff like having data
    // for the same channel in multiple packets.
   do {
      stbi__pic_packet *packet;

      if (num_packets==sizeof(packets)/sizeof(packets[0]))
         return stbi__errpuc("bad format","too many packets");

      packet = &packets[num_packets++];

      chained = stbi__get8(s);
      packet->size    = stbi__get8(s);
      packet->type    = stbi__get8(s);
      packet->channel = stbi__get8(s);

      act_comp |= packet->channel;

      if (stbi__at_eof(s))          return stbi__errpuc("bad file","file too short (reading packets)");
      if (packet->size != 8)  return stbi__errpuc("bad format","packet isn't 8bpp");
   } while (chained);

   *comp = (act_comp & 0x10 ? 4 : 3); // has alpha channel?

   for(y=0; y<height; ++y) {
      int packet_idx;

      for(packet_idx=0; packet_idx < num_packets; ++packet_idx) {
         stbi__pic_packet *packet = &packets[packet_idx];
         stbi_uc *dest = result+y*width*4;

         switch (packet->type) {
            default:
               return stbi__errpuc("bad format","packet has bad compression type");

            case 0: {//uncompressed
               int x;

               for(x=0;x<width;++x, dest+=4)
                  if (!stbi__readval(s,packet->channel,dest))
                     return 0;
               break;
            }

            case 1://Pure RLE
               {
                  int left=width, i;

                  while (left>0) {
                     stbi_uc count,value[4];

                     count=stbi__get8(s);
                     if (stbi__at_eof(s))   return stbi__errpuc("bad file","file too short (pure read count)");

                     if (count > left)
                        count = (stbi_uc) left;

                     if (!stbi__readval(s,packet->channel,value))  return 0;

                     for(i=0; i<count; ++i,dest+=4)
                        stbi__copyval(packet->channel,dest,value);
                     left -= count;
                  }
               }
               break;

            case 2: {//Mixed RLE
               int left=width;
               while (left>0) {
                  int count = stbi__get8(s), i;
                  if (stbi__at_eof(s))  return stbi__errpuc("bad file","file too short (mixed read count)");

                  if (count >= 128) { // Repeated
                     stbi_uc value[4];
                     int i;

                     if (count==128)
                        count = stbi__get16be(s);
                     else
                        count -= 127;
                     if (count > left)
                        return stbi__errpuc("bad file","scanline overrun");

                     if (!stbi__readval(s,packet->channel,value))
                        return 0;

                     for(i=0;i<count;++i, dest += 4)
                        stbi__copyval(packet->channel,dest,value);
                  } else { // Raw
                     ++count;
                     if (count>left) return stbi__errpuc("bad file","scanline overrun");

                     for(i=0;i<count;++i, dest+=4)
                        if (!stbi__readval(s,packet->channel,dest))
                           return 0;
                  }
                  left-=count;
               }
               break;
            }
         }
      }
   }

   return result;
}

static stbi_uc *stbi__pic_load(stbi__context *s,int *px,int *py,int *comp,int req_comp)
{
   stbi_uc *result;
   int i, x,y;

   for (i=0; i<92; ++i)
      stbi__get8(s);

   x = stbi__get16be(s);
   y = stbi__get16be(s);
   if (stbi__at_eof(s))  return stbi__errpuc("bad file","file too short (pic header)");
   if ((1 << 28) / x < y) return stbi__errpuc("too large", "Image too large to decode");

   stbi__get32be(s); //skip `ratio'
   stbi__get16be(s); //skip `fields'
   stbi__get16be(s); //skip `pad'

   // intermediate buffer is RGBA
   result = (stbi_uc *) stbi__malloc(x*y*4);
   memset(result, 0xff, x*y*4);

   if (!stbi__pic_load_core(s,x,y,comp, result)) {
      STBI_FREE(result);
      result=0;
   }
   *px = x;
   *py = y;
   if (req_comp == 0) req_comp = *comp;
   result=stbi__convert_format(result,4,req_comp,x,y);

   return result;
}

static int stbi__pic_test(stbi__context *s)
{
   int r = stbi__pic_test_core(s);
   stbi__rewind(s);
   return r;
}
#endif

// *************************************************************************************************
// GIF loader -- public domain by Jean-Marc Lienher -- simplified/shrunk by stb

#ifndef STBI_NO_GIF
typedef struct
{
   stbi__int16 prefix;
   stbi_uc first;
   stbi_uc suffix;
} stbi__gif_lzw;

typedef struct
{
   int w,h;
   stbi_uc *out;                 // output buffer (always 4 components)
   int flags, bgindex, ratio, transparent, eflags;
   stbi_uc  pal[256][4];
   stbi_uc lpal[256][4];
   stbi__gif_lzw codes[4096];
   stbi_uc *color_table;
   int parse, step;
   int lflags;
   int start_x, start_y;
   int max_x, max_y;
   int cur_x, cur_y;
   int line_size;
} stbi__gif;

static int stbi__gif_test_raw(stbi__context *s)
{
   int sz;
   if (stbi__get8(s) != 'G' || stbi__get8(s) != 'I' || stbi__get8(s) != 'F' || stbi__get8(s) != '8') return 0;
   sz = stbi__get8(s);
   if (sz != '9' && sz != '7') return 0;
   if (stbi__get8(s) != 'a') return 0;
   return 1;
}

static int stbi__gif_test(stbi__context *s)
{
   int r = stbi__gif_test_raw(s);
   stbi__rewind(s);
   return r;
}

static void stbi__gif_parse_colortable(stbi__context *s, stbi_uc pal[256][4], int num_entries, int transp)
{
   int i;
   for (i=0; i < num_entries; ++i) {
      pal[i][2] = stbi__get8(s);
      pal[i][1] = stbi__get8(s);
      pal[i][0] = stbi__get8(s);
      pal[i][3] = transp == i ? 0 : 255;
   }
}

static int stbi__gif_header(stbi__context *s, stbi__gif *g, int *comp, int is_info)
{
   stbi_uc version;
   if (stbi__get8(s) != 'G' || stbi__get8(s) != 'I' || stbi__get8(s) != 'F' || stbi__get8(s) != '8')
      return stbi__err("not GIF", "Corrupt GIF");

   version = stbi__get8(s);
   if (version != '7' && version != '9')    return stbi__err("not GIF", "Corrupt GIF");
   if (stbi__get8(s) != 'a')                return stbi__err("not GIF", "Corrupt GIF");

   stbi__g_failure_reason = "";
   g->w = stbi__get16le(s);
   g->h = stbi__get16le(s);
   g->flags = stbi__get8(s);
   g->bgindex = stbi__get8(s);
   g->ratio = stbi__get8(s);
   g->transparent = -1;

   if (comp != 0) *comp = 4;  // can't actually tell whether it's 3 or 4 until we parse the comments

   if (is_info) return 1;

   if (g->flags & 0x80)
      stbi__gif_parse_colortable(s,g->pal, 2 << (g->flags & 7), -1);

   return 1;
}

static int stbi__gif_info_raw(stbi__context *s, int *x, int *y, int *comp)
{
   stbi__gif g;
   if (!stbi__gif_header(s, &g, comp, 1)) {
      stbi__rewind( s );
      return 0;
   }
   if (x) *x = g.w;
   if (y) *y = g.h;
   return 1;
}

static void stbi__out_gif_code(stbi__gif *g, stbi__uint16 code)
{
   stbi_uc *p, *c;

   // recurse to decode the prefixes, since the linked-list is backwards,
   // and working backwards through an interleaved image would be nasty
   if (g->codes[code].prefix >= 0)
      stbi__out_gif_code(g, g->codes[code].prefix);

   if (g->cur_y >= g->max_y) return;

   p = &g->out[g->cur_x + g->cur_y];
   c = &g->color_table[g->codes[code].suffix * 4];

   if (c[3] >= 128) {
      p[0] = c[2];
      p[1] = c[1];
      p[2] = c[0];
      p[3] = c[3];
   }
   g->cur_x += 4;

   if (g->cur_x >= g->max_x) {
      g->cur_x = g->start_x;
      g->cur_y += g->step;

      while (g->cur_y >= g->max_y && g->parse > 0) {
         g->step = (1 << g->parse) * g->line_size;
         g->cur_y = g->start_y + (g->step >> 1);
         --g->parse;
      }
   }
}

static stbi_uc *stbi__process_gif_raster(stbi__context *s, stbi__gif *g)
{
   stbi_uc lzw_cs;
   stbi__int32 len, code;
   stbi__uint32 first;
   stbi__int32 codesize, codemask, avail, oldcode, bits, valid_bits, clear;
   stbi__gif_lzw *p;

   lzw_cs = stbi__get8(s);
   if (lzw_cs > 12) return NULL;
   clear = 1 << lzw_cs;
   first = 1;
   codesize = lzw_cs + 1;
   codemask = (1 << codesize) - 1;
   bits = 0;
   valid_bits = 0;
   for (code = 0; code < clear; code++) {
      g->codes[code].prefix = -1;
      g->codes[code].first = (stbi_uc) code;
      g->codes[code].suffix = (stbi_uc) code;
   }

   // support no starting clear code
   avail = clear+2;
   oldcode = -1;

   len = 0;
   for(;;) {
      if (valid_bits < codesize) {
         if (len == 0) {
            len = stbi__get8(s); // start new block
            if (len == 0)
               return g->out;
         }
         --len;
         bits |= (stbi__int32) stbi__get8(s) << valid_bits;
         valid_bits += 8;
      } else {
         stbi__int32 code = bits & codemask;
         bits >>= codesize;
         valid_bits -= codesize;
         // @OPTIMIZE: is there some way we can accelerate the non-clear path?
         if (code == clear) {  // clear code
            codesize = lzw_cs + 1;
            codemask = (1 << codesize) - 1;
            avail = clear + 2;
            oldcode = -1;
            first = 0;
         } else if (code == clear + 1) { // end of stream code
            stbi__skip(s, len);
            while ((len = stbi__get8(s)) > 0)
               stbi__skip(s,len);
            return g->out;
         } else if (code <= avail) {
            if (first) return stbi__errpuc("no clear code", "Corrupt GIF");

            if (oldcode >= 0) {
               p = &g->codes[avail++];
               if (avail > 4096)        return stbi__errpuc("too many codes", "Corrupt GIF");
               p->prefix = (stbi__int16) oldcode;
               p->first = g->codes[oldcode].first;
               p->suffix = (code == avail) ? p->first : g->codes[code].first;
            } else if (code == avail)
               return stbi__errpuc("illegal code in raster", "Corrupt GIF");

            stbi__out_gif_code(g, (stbi__uint16) code);

            if ((avail & codemask) == 0 && avail <= 0x0FFF) {
               codesize++;
               codemask = (1 << codesize) - 1;
            }

            oldcode = code;
         } else {
            return stbi__errpuc("illegal code in raster", "Corrupt GIF");
         }
      }
   }
}

static void stbi__fill_gif_background(stbi__gif *g)
{
   int i;
   stbi_uc *c = g->pal[g->bgindex];
   // @OPTIMIZE: write a dword at a time
   for (i = 0; i < g->w * g->h * 4; i += 4) {
      stbi_uc *p  = &g->out[i];
      p[0] = c[2];
      p[1] = c[1];
      p[2] = c[0];
      p[3] = c[3];
   }
}

// this function is designed to support animated gifs, although stb_image doesn't support it
static stbi_uc *stbi__gif_load_next(stbi__context *s, stbi__gif *g, int *comp, int req_comp)
{
   int i;
   stbi_uc *old_out = 0;

   if (g->out == 0) {
      if (!stbi__gif_header(s, g, comp,0))     return 0; // stbi__g_failure_reason set by stbi__gif_header
      g->out = (stbi_uc *) stbi__malloc(4 * g->w * g->h);
      if (g->out == 0)                      return stbi__errpuc("outofmem", "Out of memory");
      stbi__fill_gif_background(g);
   } else {
      // animated-gif-only path
      if (((g->eflags & 0x1C) >> 2) == 3) {
         old_out = g->out;
         g->out = (stbi_uc *) stbi__malloc(4 * g->w * g->h);
         if (g->out == 0)                   return stbi__errpuc("outofmem", "Out of memory");
         memcpy(g->out, old_out, g->w*g->h*4);
      }
   }

   for (;;) {
      switch (stbi__get8(s)) {
         case 0x2C: /* Image Descriptor */
         {
            stbi__int32 x, y, w, h;
            stbi_uc *o;

            x = stbi__get16le(s);
            y = stbi__get16le(s);
            w = stbi__get16le(s);
            h = stbi__get16le(s);
            if (((x + w) > (g->w)) || ((y + h) > (g->h)))
               return stbi__errpuc("bad Image Descriptor", "Corrupt GIF");

            g->line_size = g->w * 4;
            g->start_x = x * 4;
            g->start_y = y * g->line_size;
            g->max_x   = g->start_x + w * 4;
            g->max_y   = g->start_y + h * g->line_size;
            g->cur_x   = g->start_x;
            g->cur_y   = g->start_y;

            g->lflags = stbi__get8(s);

            if (g->lflags & 0x40) {
               g->step = 8 * g->line_size; // first interlaced spacing
               g->parse = 3;
            } else {
               g->step = g->line_size;
               g->parse = 0;
            }

            if (g->lflags & 0x80) {
               stbi__gif_parse_colortable(s,g->lpal, 2 << (g->lflags & 7), g->eflags & 0x01 ? g->transparent : -1);
               g->color_table = (stbi_uc *) g->lpal;
            } else if (g->flags & 0x80) {
               for (i=0; i < 256; ++i)  // @OPTIMIZE: stbi__jpeg_reset only the previous transparent
                  g->pal[i][3] = 255;
               if (g->transparent >= 0 && (g->eflags & 0x01))
                  g->pal[g->transparent][3] = 0;
               g->color_table = (stbi_uc *) g->pal;
            } else
               return stbi__errpuc("missing color table", "Corrupt GIF");

            o = stbi__process_gif_raster(s, g);
            if (o == NULL) return NULL;

            if (req_comp && req_comp != 4)
               o = stbi__convert_format(o, 4, req_comp, g->w, g->h);
            return o;
         }

         case 0x21: // Comment Extension.
         {
            int len;
            if (stbi__get8(s) == 0xF9) { // Graphic Control Extension.
               len = stbi__get8(s);
               if (len == 4) {
                  g->eflags = stbi__get8(s);
                  stbi__get16le(s); // delay
                  g->transparent = stbi__get8(s);
               } else {
                  stbi__skip(s, len);
                  break;
               }
            }
            while ((len = stbi__get8(s)) != 0)
               stbi__skip(s, len);
            break;
         }

         case 0x3B: // gif stream termination code
            return (stbi_uc *) s; // using '1' causes warning on some compilers

         default:
            return stbi__errpuc("unknown code", "Corrupt GIF");
      }
   }
}

static stbi_uc *stbi__gif_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   stbi_uc *u = 0;
   stbi__gif g;
   memset(&g, 0, sizeof(g));

   u = stbi__gif_load_next(s, &g, comp, req_comp);
   if (u == (stbi_uc *) s) u = 0;  // end of animated gif marker
   if (u) {
      *x = g.w;
      *y = g.h;
   }

   return u;
}

static int stbi__gif_info(stbi__context *s, int *x, int *y, int *comp)
{
   return stbi__gif_info_raw(s,x,y,comp);
}
#endif

// *************************************************************************************************
// Radiance RGBE HDR loader
// originally by Nicolas Schulz
#ifndef STBI_NO_HDR
static int stbi__hdr_test_core(stbi__context *s)
{
   const char *signature = "#?RADIANCE\n";
   int i;
   for (i=0; signature[i]; ++i)
      if (stbi__get8(s) != signature[i])
         return 0;
   return 1;
}

static int stbi__hdr_test(stbi__context* s)
{
   int r = stbi__hdr_test_core(s);
   stbi__rewind(s);
   return r;
}

#define STBI__HDR_BUFLEN  1024
static char *stbi__hdr_gettoken(stbi__context *z, char *buffer)
{
   int len=0;
   char c = '\0';

   c = (char) stbi__get8(z);

   while (!stbi__at_eof(z) && c != '\n') {
      buffer[len++] = c;
      if (len == STBI__HDR_BUFLEN-1) {
         // flush to end of line
         while (!stbi__at_eof(z) && stbi__get8(z) != '\n')
            ;
         break;
      }
      c = (char) stbi__get8(z);
   }

   buffer[len] = 0;
   return buffer;
}

static void stbi__hdr_convert(float *output, stbi_uc *input, int req_comp)
{
   if ( input[3] != 0 ) {
      float f1;
      // Exponent
      f1 = (float) ldexp(1.0f, input[3] - (int)(128 + 8));
      if (req_comp <= 2)
         output[0] = (input[0] + input[1] + input[2]) * f1 / 3;
      else {
         output[0] = input[0] * f1;
         output[1] = input[1] * f1;
         output[2] = input[2] * f1;
      }
      if (req_comp == 2) output[1] = 1;
      if (req_comp == 4) output[3] = 1;
   } else {
      switch (req_comp) {
         case 4: output[3] = 1; /* fallthrough */
         case 3: output[0] = output[1] = output[2] = 0;
                 break;
         case 2: output[1] = 1; /* fallthrough */
         case 1: output[0] = 0;
                 break;
      }
   }
}

static float *stbi__hdr_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   char buffer[STBI__HDR_BUFLEN];
   char *token;
   int valid = 0;
   int width, height;
   stbi_uc *scanline;
   float *hdr_data;
   int len;
   unsigned char count, value;
   int i, j, k, c1,c2, z;


   // Check identifier
   if (strcmp(stbi__hdr_gettoken(s,buffer), "#?RADIANCE") != 0)
      return stbi__errpf("not HDR", "Corrupt HDR image");

   // Parse header
   for(;;) {
      token = stbi__hdr_gettoken(s,buffer);
      if (token[0] == 0) break;
      if (strcmp(token, "FORMAT=32-bit_rle_rgbe") == 0) valid = 1;
   }

   if (!valid)    return stbi__errpf("unsupported format", "Unsupported HDR format");

   // Parse width and height
   // can't use sscanf() if we're not using stdio!
   token = stbi__hdr_gettoken(s,buffer);
   if (strncmp(token, "-Y ", 3))  return stbi__errpf("unsupported data layout", "Unsupported HDR format");
   token += 3;
   height = (int) strtol(token, &token, 10);
   while (*token == ' ') ++token;
   if (strncmp(token, "+X ", 3))  return stbi__errpf("unsupported data layout", "Unsupported HDR format");
   token += 3;
   width = (int) strtol(token, NULL, 10);

   *x = width;
   *y = height;

   if (comp) *comp = 3;
   if (req_comp == 0) req_comp = 3;

   // Read data
   hdr_data = (float *) stbi__malloc(height * width * req_comp * sizeof(float));

   // Load image data
   // image data is stored as some number of sca
   if ( width < 8 || width >= 32768) {
      // Read flat data
      for (j=0; j < height; ++j) {
         for (i=0; i < width; ++i) {
            stbi_uc rgbe[4];
           main_decode_loop:
            stbi__getn(s, rgbe, 4);
            stbi__hdr_convert(hdr_data + j * width * req_comp + i * req_comp, rgbe, req_comp);
         }
      }
   } else {
      // Read RLE-encoded data
      scanline = NULL;

      for (j = 0; j < height; ++j) {
         c1 = stbi__get8(s);
         c2 = stbi__get8(s);
         len = stbi__get8(s);
         if (c1 != 2 || c2 != 2 || (len & 0x80)) {
            // not run-length encoded, so we have to actually use THIS data as a decoded
            // pixel (note this can't be a valid pixel--one of RGB must be >= 128)
            stbi_uc rgbe[4];
            rgbe[0] = (stbi_uc) c1;
            rgbe[1] = (stbi_uc) c2;
            rgbe[2] = (stbi_uc) len;
            rgbe[3] = (stbi_uc) stbi__get8(s);
            stbi__hdr_convert(hdr_data, rgbe, req_comp);
            i = 1;
            j = 0;
            STBI_FREE(scanline);
            goto main_decode_loop; // yes, this makes no sense
         }
         len <<= 8;
         len |= stbi__get8(s);
         if (len != width) { STBI_FREE(hdr_data); STBI_FREE(scanline); return stbi__errpf("invalid decoded scanline length", "corrupt HDR"); }
         if (scanline == NULL) scanline = (stbi_uc *) stbi__malloc(width * 4);

         for (k = 0; k < 4; ++k) {
            i = 0;
            while (i < width) {
               count = stbi__get8(s);
               if (count > 128) {
                  // Run
                  value = stbi__get8(s);
                  count -= 128;
                  for (z = 0; z < count; ++z)
                     scanline[i++ * 4 + k] = value;
               } else {
                  // Dump
                  for (z = 0; z < count; ++z)
                     scanline[i++ * 4 + k] = stbi__get8(s);
               }
            }
         }
         for (i=0; i < width; ++i)
            stbi__hdr_convert(hdr_data+(j*width + i)*req_comp, scanline + i*4, req_comp);
      }
      STBI_FREE(scanline);
   }

   return hdr_data;
}

static int stbi__hdr_info(stbi__context *s, int *x, int *y, int *comp)
{
   char buffer[STBI__HDR_BUFLEN];
   char *token;
   int valid = 0;

   if (strcmp(stbi__hdr_gettoken(s,buffer), "#?RADIANCE") != 0) {
       stbi__rewind( s );
       return 0;
   }

   for(;;) {
      token = stbi__hdr_gettoken(s,buffer);
      if (token[0] == 0) break;
      if (strcmp(token, "FORMAT=32-bit_rle_rgbe") == 0) valid = 1;
   }

   if (!valid) {
       stbi__rewind( s );
       return 0;
   }
   token = stbi__hdr_gettoken(s,buffer);
   if (strncmp(token, "-Y ", 3)) {
       stbi__rewind( s );
       return 0;
   }
   token += 3;
   *y = (int) strtol(token, &token, 10);
   while (*token == ' ') ++token;
   if (strncmp(token, "+X ", 3)) {
       stbi__rewind( s );
       return 0;
   }
   token += 3;
   *x = (int) strtol(token, NULL, 10);
   *comp = 3;
   return 1;
}
#endif // STBI_NO_HDR

#ifndef STBI_NO_BMP
static int stbi__bmp_info(stbi__context *s, int *x, int *y, int *comp)
{
   int hsz;
   if (stbi__get8(s) != 'B' || stbi__get8(s) != 'M') {
       stbi__rewind( s );
       return 0;
   }
   stbi__skip(s,12);
   hsz = stbi__get32le(s);
   if (hsz != 12 && hsz != 40 && hsz != 56 && hsz != 108 && hsz != 124) {
       stbi__rewind( s );
       return 0;
   }
   if (hsz == 12) {
      *x = stbi__get16le(s);
      *y = stbi__get16le(s);
   } else {
      *x = stbi__get32le(s);
      *y = stbi__get32le(s);
   }
   if (stbi__get16le(s) != 1) {
       stbi__rewind( s );
       return 0;
   }
   *comp = stbi__get16le(s) / 8;
   return 1;
}
#endif

#ifndef STBI_NO_PSD
static int stbi__psd_info(stbi__context *s, int *x, int *y, int *comp)
{
   int channelCount;
   if (stbi__get32be(s) != 0x38425053) {
       stbi__rewind( s );
       return 0;
   }
   if (stbi__get16be(s) != 1) {
       stbi__rewind( s );
       return 0;
   }
   stbi__skip(s, 6);
   channelCount = stbi__get16be(s);
   if (channelCount < 0 || channelCount > 16) {
       stbi__rewind( s );
       return 0;
   }
   *y = stbi__get32be(s);
   *x = stbi__get32be(s);
   if (stbi__get16be(s) != 8) {
       stbi__rewind( s );
       return 0;
   }
   if (stbi__get16be(s) != 3) {
       stbi__rewind( s );
       return 0;
   }
   *comp = 4;
   return 1;
}
#endif

#ifndef STBI_NO_PIC
static int stbi__pic_info(stbi__context *s, int *x, int *y, int *comp)
{
   int act_comp=0,num_packets=0,chained;
   stbi__pic_packet packets[10];

   stbi__skip(s, 92);

   *x = stbi__get16be(s);
   *y = stbi__get16be(s);
   if (stbi__at_eof(s))  return 0;
   if ( (*x) != 0 && (1 << 28) / (*x) < (*y)) {
       stbi__rewind( s );
       return 0;
   }

   stbi__skip(s, 8);

   do {
      stbi__pic_packet *packet;

      if (num_packets==sizeof(packets)/sizeof(packets[0]))
         return 0;

      packet = &packets[num_packets++];
      chained = stbi__get8(s);
      packet->size    = stbi__get8(s);
      packet->type    = stbi__get8(s);
      packet->channel = stbi__get8(s);
      act_comp |= packet->channel;

      if (stbi__at_eof(s)) {
          stbi__rewind( s );
          return 0;
      }
      if (packet->size != 8) {
          stbi__rewind( s );
          return 0;
      }
   } while (chained);

   *comp = (act_comp & 0x10 ? 4 : 3);

   return 1;
}
#endif

// *************************************************************************************************
// Portable Gray Map and Portable Pixel Map loader
// by Ken Miller
//
// PGM: http://netpbm.sourceforge.net/doc/pgm.html
// PPM: http://netpbm.sourceforge.net/doc/ppm.html
//
// Known limitations:
//    Does not support comments in the header section
//    Does not support ASCII image data (formats P2 and P3)
//    Does not support 16-bit-per-channel

#ifndef STBI_NO_PNM

static int      stbi__pnm_test(stbi__context *s)
{
   char p, t;
   p = (char) stbi__get8(s);
   t = (char) stbi__get8(s);
   if (p != 'P' || (t != '5' && t != '6')) {
       stbi__rewind( s );
       return 0;
   }
   return 1;
}

static stbi_uc *stbi__pnm_load(stbi__context *s, int *x, int *y, int *comp, int req_comp)
{
   stbi_uc *out;
   if (!stbi__pnm_info(s, (int *)&s->img_x, (int *)&s->img_y, (int *)&s->img_n))
      return 0;
   *x = s->img_x;
   *y = s->img_y;
   *comp = s->img_n;

   out = (stbi_uc *) stbi__malloc(s->img_n * s->img_x * s->img_y);
   if (!out) return stbi__errpuc("outofmem", "Out of memory");
   stbi__getn(s, out, s->img_n * s->img_x * s->img_y);

   if (req_comp && req_comp != s->img_n) {
      out = stbi__convert_format(out, s->img_n, req_comp, s->img_x, s->img_y);
      if (out == NULL) return out; // stbi__convert_format frees input on failure
   }
   return out;
}

static int      stbi__pnm_isspace(char c)
{
   return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}

static void     stbi__pnm_skip_whitespace(stbi__context *s, char *c)
{
   while (!stbi__at_eof(s) && stbi__pnm_isspace(*c))
      *c = (char) stbi__get8(s);
}

static int      stbi__pnm_isdigit(char c)
{
   return c >= '0' && c <= '9';
}

static int      stbi__pnm_getinteger(stbi__context *s, char *c)
{
   int value = 0;

   while (!stbi__at_eof(s) && stbi__pnm_isdigit(*c)) {
      value = value*10 + (*c - '0');
      *c = (char) stbi__get8(s);
   }

   return value;
}

static int      stbi__pnm_info(stbi__context *s, int *x, int *y, int *comp)
{
   int maxv;
   char c, p, t;

   stbi__rewind( s );

   // Get identifier
   p = (char) stbi__get8(s);
   t = (char) stbi__get8(s);
   if (p != 'P' || (t != '5' && t != '6')) {
       stbi__rewind( s );
       return 0;
   }

   *comp = (t == '6') ? 3 : 1;  // '5' is 1-component .pgm; '6' is 3-component .ppm

   c = (char) stbi__get8(s);
   stbi__pnm_skip_whitespace(s, &c);

   *x = stbi__pnm_getinteger(s, &c); // read width
   stbi__pnm_skip_whitespace(s, &c);

   *y = stbi__pnm_getinteger(s, &c); // read height
   stbi__pnm_skip_whitespace(s, &c);

   maxv = stbi__pnm_getinteger(s, &c);  // read max value

   if (maxv > 255)
      return stbi__err("max value > 255", "PPM image not 8-bit");
   else
      return 1;
}
#endif

static int stbi__info_main(stbi__context *s, int *x, int *y, int *comp)
{
   #ifndef STBI_NO_JPEG
   if (stbi__jpeg_info(s, x, y, comp)) return 1;
   #endif

   #ifndef STBI_NO_PNG
   if (stbi__png_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_GIF
   if (stbi__gif_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_BMP
   if (stbi__bmp_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_PSD
   if (stbi__psd_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_PIC
   if (stbi__pic_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_PNM
   if (stbi__pnm_info(s, x, y, comp))  return 1;
   #endif

   #ifndef STBI_NO_HDR
   if (stbi__hdr_info(s, x, y, comp))  return 1;
   #endif

   // test tga last because it's a crappy test!
   #ifndef STBI_NO_TGA
   if (stbi__tga_info(s, x, y, comp))
       return 1;
   #endif
   return stbi__err("unknown image type", "Image not of any known type, or corrupt");
}

#ifndef STBI_NO_STDIO
STBIDEF int stbi_info(char *filename, int *x, int *y, int *comp)
{
    FILE *f = stbi__fopen(filename, "rb");
    int result;
    if (!f) return stbi__err("can't fopen", "Unable to open file");
    result = stbi_info_from_file(f, x, y, comp);
    fclose(f);
    return result;
}

STBIDEF int stbi_info_from_file(FILE *f, int *x, int *y, int *comp)
{
   int r;
   stbi__context s;
   long pos = ftell(f);
   stbi__start_file(&s, f);
   r = stbi__info_main(&s,x,y,comp);
   fseek(f,pos,SEEK_SET);
   return r;
}
#endif // !STBI_NO_STDIO

STBIDEF int stbi_info_from_memory(stbi_uc *buffer, int len, int *x, int *y, int *comp)
{
   stbi__context s;
   stbi__start_mem(&s,buffer,len);
   return stbi__info_main(&s,x,y,comp);
}

STBIDEF int stbi_info_from_callbacks(stbi_io_callbacks *c, void *user, int *x, int *y, int *comp)
{
   stbi__context s;
   stbi__start_callbacks(&s, (stbi_io_callbacks *) c, user);
   return stbi__info_main(&s,x,y,comp);
}

/*
   revision history:
      2.06  (2015-04-19) fix bug where PSD returns wrong '*comp' value
      2.05  (2015-04-19) fix bug in progressive JPEG handling, fix warning
      2.04  (2015-04-15) try to re-enable SIMD on MinGW 64-bit
      2.03  (2015-04-12) extra corruption checking (mmozeiko)
                         stbi_set_flip_vertically_on_load (nguillemot)
                         fix NEON support; fix mingw support
      2.02  (2015-01-19) fix incorrect assert, fix warning
      2.01  (2015-01-17) fix various warnings; suppress SIMD on gcc 32-bit without -msse2
      2.00b (2014-12-25) fix STBI_MALLOC in progressive JPEG
      2.00  (2014-12-25) optimize JPG, including x86 SSE2 & NEON SIMD (ryg)
                         progressive JPEG (stb)
                         PGM/PPM support (Ken Miller)
                         STBI_MALLOC,STBI_REALLOC,STBI_FREE
                         GIF bugfix -- seemingly never worked
                         STBI_NO_*, STBI_ONLY_*
      1.48  (2014-12-14) fix incorrectly-named assert()
      1.47  (2014-12-14) 1/2/4-bit PNG support, both direct and paletted (Omar Cornut & stb)
                         optimize PNG (ryg)
                         fix bug in interlaced PNG with user-specified channel count (stb)
      1.46  (2014-08-26)
              fix broken tRNS chunk (colorkey-style transparency) in non-paletted PNG
      1.45  (2014-08-16)
              fix MSVC-ARM internal compiler error by wrapping malloc
      1.44  (2014-08-07)
              various warning fixes from Ronny Chevalier
      1.43  (2014-07-15)
              fix MSVC-only compiler problem in code changed in 1.42
      1.42  (2014-07-09)
              don't define _CRT_SECURE_NO_WARNINGS (affects user code)
              fixes to stbi__cleanup_jpeg path
              added STBI_ASSERT to avoid requiring assert.h
      1.41  (2014-06-25)
              fix search&replace from 1.36 that messed up comments/error messages
      1.40  (2014-06-22)
              fix gcc struct-initialization warning
      1.39  (2014-06-15)
              fix to TGA optimization when req_comp != number of components in TGA;
              fix to GIF loading because BMP wasn't rewinding (whoops, no GIFs in my test suite)
              add support for BMP version 5 (more ignored fields)
      1.38  (2014-06-06)
              suppress MSVC warnings on integer casts truncating values
              fix accidental rename of 'skip' field of I/O
      1.37  (2014-06-04)
              remove duplicate typedef
      1.36  (2014-06-03)
              convert to header file single-file library
              if de-iphone isn't set, load iphone images color-swapped instead of returning NULL
      1.35  (2014-05-27)
              various warnings
              fix broken STBI_SIMD path
              fix bug where stbi_load_from_file no longer left file pointer in correct place
              fix broken non-easy path for 32-bit BMP (possibly never used)
              TGA optimization by Arseny Kapoulkine
      1.34  (unknown)
              use STBI_NOTUSED in stbi__resample_row_generic(), fix one more leak in tga failure case
      1.33  (2011-07-14)
              make stbi_is_hdr work in STBI_NO_HDR (as specified), minor compiler-friendly improvements
      1.32  (2011-07-13)
              support for "info" function for all supported filetypes (SpartanJ)
      1.31  (2011-06-20)
              a few more leak fixes, bug in PNG handling (SpartanJ)
      1.30  (2011-06-11)
              added ability to load files via callbacks to accomidate custom input streams (Ben Wenger)
              removed deprecated format-specific test/load functions
              removed support for installable file formats (stbi_loader) -- would have been broken for IO callbacks anyway
              error cases in bmp and tga give messages and don't leak (Raymond Barbiero, grisha)
              fix inefficiency in decoding 32-bit BMP (David Woo)
      1.29  (2010-08-16)
              various warning fixes from Aurelien Pocheville
      1.28  (2010-08-01)
              fix bug in GIF palette transparency (SpartanJ)
      1.27  (2010-08-01)
              cast-to-stbi_uc to fix warnings
      1.26  (2010-07-24)
              fix bug in file buffering for PNG reported by SpartanJ
      1.25  (2010-07-17)
              refix trans_data warning (Won Chun)
      1.24  (2010-07-12)
              perf improvements reading from files on platforms with lock-heavy fgetc()
              minor perf improvements for jpeg
              deprecated type-specific functions so we'll get feedback if they're needed
              attempt to fix trans_data warning (Won Chun)
      1.23    fixed bug in iPhone support
      1.22  (2010-07-10)
              removed image *writing* support
              stbi_info support from Jetro Lauha
              GIF support from Jean-Marc Lienher
              iPhone PNG-extensions from James Brown
              warning-fixes from Nicolas Schulz and Janez Zemva (i.stbi__err. Janez (U+017D)emva)
      1.21    fix use of 'stbi_uc' in header (reported by jon blow)
      1.20    added support for Softimage PIC, by Tom Seddon
      1.19    bug in interlaced PNG corruption check (found by ryg)
      1.18  (2008-08-02)
              fix a threading bug (local mutable static)
      1.17    support interlaced PNG
      1.16    major bugfix - stbi__convert_format converted one too many pixels
      1.15    initialize some fields for thread safety
      1.14    fix threadsafe conversion bug
              header-file-only version (#define STBI_HEADER_FILE_ONLY before including)
      1.13    threadsafe
      1.12    const qualifiers in the API
      1.11    Support installable IDCT, colorspace conversion routines
      1.10    Fixes for 64-bit (don't use "unsigned long")
              optimized upsampling by Fabian "ryg" Giesen
      1.09    Fix format-conversion for PSD code (bad global variables!)
      1.08    Thatcher Ulrich's PSD code integrated by Nicolas Schulz
      1.07    attempt to fix C++ warning/errors again
      1.06    attempt to fix C++ warning/errors again
      1.05    fix TGA loading to return correct *comp and use good luminance calc
      1.04    default float alpha is 1, not 255; use 'void *' for stbi_image_free
      1.03    bugfixes to STBI_NO_STDIO, STBI_NO_HDR
      1.02    support for (subset of) HDR files, float interface for preferred access to them
      1.01    fix bug: possible bug in handling right-side up bmps... not sure
              fix bug: the stbi__bmp_load() and stbi__tga_load() functions didn't work at all
      1.00    interface to zlib that skips zlib header
      0.99    correct handling of alpha in palette
      0.98    TGA loader by lonesock; dynamically add loaders (untested)
      0.97    jpeg errors on too large a file; also catch another malloc failure
      0.96    fix detection of invalid v value - particleman@mollyrocket forum
      0.95    during header scan, seek to markers in case of padding
      0.94    STBI_NO_STDIO to disable stdio usage; rename all #defines the same
      0.93    handle jpegtran output; verbose errors
      0.92    read 4,8,16,24,32-bit BMP files of several formats
      0.91    output 24-bit Windows 3.0 BMP files
      0.90    fix a few more warnings; bump version number to approach 1.0
      0.61    bugfixes due to Marc LeBlanc, Christopher Lloyd
      0.60    fix compiling as c++
      0.59    fix warnings: merge Dave Moore's -Wall fixes
      0.58    fix bug: zlib uncompressed mode len/nlen was wrong endian
      0.57    fix bug: jpg last huffman symbol before marker was >9 bits but less than 16 available
      0.56    fix bug: zlib uncompressed mode len vs. nlen
      0.55    fix bug: restart_interval not initialized to 0
      0.54    allow NULL for 'int *comp'
      0.53    fix bug in png 3->4; speedup png decoding
      0.52    png handles req_comp=3,4 directly; minor cleanup; jpeg comments
      0.51    obey req_comp requests, 1-component jpegs return as 1-component,
              on 'test' only check type, not whether we support this variant
      0.50  (2006-11-19)
              first released version
*/

/* stb_image_write - v0.98 - public domain - http://nothings.org/stb/stb_image_write.h
   writes out PNG/BMP/TGA images to C stdio - Sean Barrett 2010
                            no warranty implied; use at your own risk


   Before #including,

       #define STB_IMAGE_WRITE_IMPLEMENTATION

   in the file that you want to have the implementation.

   Will probably not work correctly with strict-aliasing optimizations.

ABOUT:

   This header file is a library for writing images to C stdio. It could be
   adapted to write to memory or a general streaming interface; let me know.

   The PNG output is not optimal; it is 20-50% larger than the file
   written by a decent optimizing implementation. This library is designed
   for source code compactness and simplicitly, not optimal image file size
   or run-time performance.

BUILDING:

   You can #define STBIW_ASSERT(x) before the #include to avoid using assert.h.
   You can #define STBIW_MALLOC(), STBIW_REALLOC(), and STBIW_FREE() to replace
   malloc,realloc,free.
   You can define STBIW_MEMMOVE() to replace memmove()

USAGE:

   There are four functions, one for each image file format:

     int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
     int stbi_write_bmp(char const *filename, int w, int h, int comp, const void *data);
     int stbi_write_tga(char const *filename, int w, int h, int comp, const void *data);
     int stbi_write_hdr(char const *filename, int w, int h, int comp, const void *data);

   Each function returns 0 on failure and non-0 on success.

   The functions create an image file defined by the parameters. The image
   is a rectangle of pixels stored from left-to-right, top-to-bottom.
   Each pixel contains 'comp' channels of data stored interleaved with 8-bits
   per channel, in the following order: 1=Y, 2=YA, 3=RGB, 4=RGBA. (Y is
   monochrome color.) The rectangle is 'w' pixels wide and 'h' pixels tall.
   The *data pointer points to the first byte of the top-left-most pixel.
   For PNG, "stride_in_bytes" is the distance in bytes from the first byte of
   a row of pixels to the first byte of the next row of pixels.

   PNG creates output files with the same number of components as the input.
   The BMP format expands Y to RGB in the file format and does not
   output alpha.

   PNG supports writing rectangles of data even when the bytes storing rows of
   data are not consecutive in memory (e.g. sub-rectangles of a larger image),
   by supplying the stride between the beginning of adjacent rows. The other
   formats do not. (Thus you cannot write a native-format BMP through the BMP
   writer, both because it is in BGR order and because it may have padding
   at the end of the line.)

   HDR expects linear float data. Since the format is always 32-bit rgb(e)
   data, alpha (if provided) is discarded, and for monochrome data it is
   replicated across all three channels.

CREDITS:

   PNG/BMP/TGA
      Sean Barrett
   HDR
      Baldur Karlsson
   TGA monochrome:
      Jean-Sebastien Guay
   misc enhancements:
      Tim Kelsey
   bugfixes:
      github:Chribba
*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#if defined(STBIW_MALLOC) && defined(STBIW_FREE) && defined(STBIW_REALLOC)
// ok
#elif !defined(STBIW_MALLOC) && !defined(STBIW_FREE) && !defined(STBIW_REALLOC)
// ok
#else
#error "Must define all or none of STBIW_MALLOC, STBIW_FREE, and STBIW_REALLOC."
#endif

#ifndef STBIW_MALLOC
#define STBIW_MALLOC(sz)    malloc(sz)
#define STBIW_REALLOC(p,sz) realloc(p,sz)
#define STBIW_FREE(p)       free(p)
#endif
#ifndef STBIW_MEMMOVE
#define STBIW_MEMMOVE(a,b,sz) memmove(a,b,sz)
#endif


#ifndef STBIW_ASSERT
#include <assert.h>
#define STBIW_ASSERT(x) assert(x)
#endif

typedef unsigned int stbiw_uint32;
typedef int stb_image_write_test[sizeof(stbiw_uint32)==4 ? 1 : -1];

static void writefv(FILE *f, const char *fmt, va_list v)
{
   while (*fmt) {
      switch (*fmt++) {
         case ' ': break;
         case '1': { unsigned char x = (unsigned char) va_arg(v, int); fputc(x,f); break; }
         case '2': { int x = va_arg(v,int); unsigned char b[2];
                     b[0] = (unsigned char) x; b[1] = (unsigned char) (x>>8);
                     fwrite(b,2,1,f); break; }
         case '4': { stbiw_uint32 x = va_arg(v,int); unsigned char b[4];
                     b[0]=(unsigned char)x; b[1]=(unsigned char)(x>>8);
                     b[2]=(unsigned char)(x>>16); b[3]=(unsigned char)(x>>24);
                     fwrite(b,4,1,f); break; }
         default:
            STBIW_ASSERT(0);
            return;
      }
   }
}

static void write3(FILE *f, unsigned char a, unsigned char b, unsigned char c)
{
   unsigned char arr[3];
   arr[0] = a, arr[1] = b, arr[2] = c;
   fwrite(arr, 3, 1, f);
}

static void write_pixels(FILE *f, int rgb_dir, int vdir, int x, int y, int comp, const void *data, int write_alpha, int scanline_pad, int expand_mono)
{
   unsigned char bg[3] = { 255, 0, 255}, px[3];
   stbiw_uint32 zero = 0;
   int i,j,k, j_end;

   if (y <= 0)
      return;

   if (vdir < 0)
      j_end = -1, j = y-1;
   else
      j_end =  y, j = 0;

   for (; j != j_end; j += vdir) {
      for (i=0; i < x; ++i) {
         const unsigned char *d = (const unsigned char *) data + (j*x+i)*comp;
         if (write_alpha < 0)
            fwrite(&d[comp-1], 1, 1, f);
         switch (comp) {
            case 1: fwrite(d, 1, 1, f);
                    break;
            case 2: if (expand_mono)
                       write3(f, d[0],d[0],d[0]); // monochrome bmp
                    else
                       fwrite(d, 1, 1, f);  // monochrome TGA
                    break;
            case 4:
               if (!write_alpha) {
                  // composite against pink background
                  for (k=0; k < 3; ++k)
                     px[k] = bg[k] + ((d[k] - bg[k]) * d[3])/255;
                  write3(f, px[1-rgb_dir],px[1],px[1+rgb_dir]);
                  break;
               }
               /* FALLTHROUGH */
            case 3:
               write3(f, d[1-rgb_dir],d[1],d[1+rgb_dir]);
               break;
         }
         if (write_alpha > 0)
            fwrite(&d[comp-1], 1, 1, f);
      }
      fwrite(&zero,scanline_pad,1,f);
   }
}

static int outfile(char const *filename, int rgb_dir, int vdir, int x, int y, int comp, int expand_mono, const void *data, int alpha, int pad, const char *fmt, ...)
{
   FILE *f;
   if (y < 0 || x < 0) return 0;
   f = fopen(filename, "wb");
   if (f) {
      va_list v;
      va_start(v, fmt);
      writefv(f, fmt, v);
      va_end(v);
      write_pixels(f,rgb_dir,vdir,x,y,comp,data,alpha,pad,expand_mono);
      fclose(f);
   }
   return f != NULL;
}

int stbi_write_bmp(char const *filename, int x, int y, int comp, const void *data)
{
   int pad = (-x*3) & 3;
   return outfile(filename,-1,-1,x,y,comp,1,data,0,pad,
           "11 4 22 4" "4 44 22 444444",
           'B', 'M', 14+40+(x*3+pad)*y, 0,0, 14+40,  // file header
            40, x,y, 1,24, 0,0,0,0,0,0);             // bitmap header
}

int stbi_write_tga(char const *filename, int x, int y, int const comp, void *data)
{
   int has_alpha = (comp == 2 || comp == 4);
   int colorbytes = has_alpha ? comp-1 : comp;
   int format = colorbytes < 2 ? 3 : 2; // 3 color channels (RGB/RGBA) = 2, 1 color channel (Y/YA) = 3
   return outfile(filename, -1,-1, x, y, comp, 0, data, has_alpha, 0,
                  "111 221 2222 11", 0,0,format, 0,0,0, 0,0,x,y, (colorbytes+has_alpha)*8, has_alpha*8);
}

// *************************************************************************************************
// Radiance RGBE HDR writer
// by Baldur Karlsson
#define stbiw__max(a, b)  ((a) > (b) ? (a) : (b))

void stbiw__linear_to_rgbe(unsigned char *rgbe, float *linear)
{
   int exponent;
   float maxcomp = stbiw__max(linear[0], stbiw__max(linear[1], linear[2]));

   if (maxcomp < 1e-32) {
      rgbe[0] = rgbe[1] = rgbe[2] = rgbe[3] = 0;
   } else {
      float normalize = (float) frexp(maxcomp, &exponent) * 256.0f/maxcomp;

      rgbe[0] = (unsigned char)(linear[0] * normalize);
      rgbe[1] = (unsigned char)(linear[1] * normalize);
      rgbe[2] = (unsigned char)(linear[2] * normalize);
      rgbe[3] = (unsigned char)(exponent + 128);
   }
}

void stbiw__write_run_data(FILE *f, int length, unsigned char databyte)
{
   unsigned char lengthbyte = (unsigned char) (length+128);
   STBIW_ASSERT(length+128 <= 255);
   fwrite(&lengthbyte, 1, 1, f);
   fwrite(&databyte, 1, 1, f);
}

void stbiw__write_dump_data(FILE *f, int length, unsigned char *data)
{
   unsigned char lengthbyte = (unsigned char )(length & 0xff);
   STBIW_ASSERT(length <= 128); // inconsistent with spec but consistent with official code
   fwrite(&lengthbyte, 1, 1, f);
   fwrite(data, length, 1, f);
}

void stbiw__write_hdr_scanline(FILE *f, int width, int comp, unsigned char *scratch, const float *scanline)
{
   unsigned char scanlineheader[4] = { 2, 2, 0, 0 };
   unsigned char rgbe[4];
   float linear[3];
   int x;

   scanlineheader[2] = (width&0xff00)>>8;
   scanlineheader[3] = (width&0x00ff);

   /* skip RLE for images too small or large */
   if (width < 8 || width >= 32768) {
      for (x=0; x < width; x++) {
         switch (comp) {
            case 4: /* fallthrough */
            case 3: linear[2] = scanline[x*comp + 2];
                    linear[1] = scanline[x*comp + 1];
                    linear[0] = scanline[x*comp + 0];
                    break;
            case 2: /* fallthrough */
            case 1: linear[0] = linear[1] = linear[2] = scanline[x*comp + 0];
                    break;
         }
         stbiw__linear_to_rgbe(rgbe, linear);
         fwrite(rgbe, 4, 1, f);
      }
   } else {
      int c,r;
      /* encode into scratch buffer */
      for (x=0; x < width; x++) {
         switch(comp) {
            case 4: /* fallthrough */
            case 3: linear[2] = scanline[x*comp + 2];
                    linear[1] = scanline[x*comp + 1];
                    linear[0] = scanline[x*comp + 0];
                    break;
            case 2: /* fallthrough */
            case 1: linear[0] = linear[1] = linear[2] = scanline[x*comp + 0];
                    break;
         }
         stbiw__linear_to_rgbe(rgbe, linear);
         scratch[x + width*0] = rgbe[0];
         scratch[x + width*1] = rgbe[1];
         scratch[x + width*2] = rgbe[2];
         scratch[x + width*3] = rgbe[3];
      }

      fwrite(scanlineheader, 4, 1, f);

      /* RLE each component separately */
      for (c=0; c < 4; c++) {
         unsigned char *comp = &scratch[width*c];

         x = 0;
         while (x < width) {
            // find first run
            r = x;
            while (r+2 < width) {
               if (comp[r] == comp[r+1] && comp[r] == comp[r+2])
                  break;
               ++r;
            }
            if (r+2 >= width)
               r = width;
            // dump up to first run
            while (x < r) {
               int len = r-x;
               if (len > 128) len = 128;
               stbiw__write_dump_data(f, len, &comp[x]);
               x += len;
            }
            // if there's a run, output it
            if (r+2 < width) { // same test as what we break out of in search loop, so only true if we break'd
               // find next byte after run
               while (r < width && comp[r] == comp[x])
                  ++r;
               // output run up to r
               while (x < r) {
                  int len = r-x;
                  if (len > 127) len = 127;
                  stbiw__write_run_data(f, len, comp[x]);
                  x += len;
               }
            }
         }
      }
   }
}

int stbi_write_hdr(char const *filename, int x, int y, int comp, const float *data)
{
   int i;
   FILE *f;
   if (y <= 0 || x <= 0 || data == NULL) return 0;
   f = fopen(filename, "wb");
   if (f) {
      /* Each component is stored separately. Allocate scratch space for full output scanline. */
      unsigned char *scratch = (unsigned char *) STBIW_MALLOC(x*4);
      fprintf(f, "#?RADIANCE\n# Written by stb_image_write.h\nFORMAT=32-bit_rle_rgbe\n"      );
      fprintf(f, "EXPOSURE=          1.0000000000000\n\n-Y %d +X %d\n"                 , y, x);
      for(i=0; i < y; i++)
         stbiw__write_hdr_scanline(f, x, comp, scratch, data + comp*i*x);
      STBIW_FREE(scratch);
      fclose(f);
   }
   return f != NULL;
}

/////////////////////////////////////////////////////////
// PNG

// stretchy buffer; stbiw__sbpush() == vector<>::push_back() -- stbiw__sbcount() == vector<>::size()
#define stbiw__sbraw(a) ((int *) (a) - 2)
#define stbiw__sbm(a)   stbiw__sbraw(a)[0]
#define stbiw__sbn(a)   stbiw__sbraw(a)[1]

#define stbiw__sbneedgrow(a,n)  ((a)==0 || stbiw__sbn(a)+n >= stbiw__sbm(a))
#define stbiw__sbmaybegrow(a,n) (stbiw__sbneedgrow(a,(n)) ? stbiw__sbgrow(a,n) : 0)
#define stbiw__sbgrow(a,n)  stbiw__sbgrowf((void **) &(a), (n), sizeof(*(a)))

#define stbiw__sbpush(a, v)      (stbiw__sbmaybegrow(a,1), (a)[stbiw__sbn(a)++] = (v))
#define stbiw__sbcount(a)        ((a) ? stbiw__sbn(a) : 0)
#define stbiw__sbfree(a)         ((a) ? STBIW_FREE(stbiw__sbraw(a)),0 : 0)

static void *stbiw__sbgrowf(void **arr, int increment, int itemsize)
{
   int m = *arr ? 2*stbiw__sbm(*arr)+increment : increment+1;
   void *p = STBIW_REALLOC(*arr ? stbiw__sbraw(*arr) : 0, itemsize * m + sizeof(int)*2);
   STBIW_ASSERT(p);
   if (p) {
      if (!*arr) ((int *) p)[1] = 0;
      *arr = (void *) ((int *) p + 2);
      stbiw__sbm(*arr) = m;
   }
   return *arr;
}

static unsigned char *stbiw__zlib_flushf(unsigned char *data, unsigned int *bitbuffer, int *bitcount)
{
   while (*bitcount >= 8) {
      stbiw__sbpush(data, (unsigned char) *bitbuffer);
      *bitbuffer >>= 8;
      *bitcount -= 8;
   }
   return data;
}

static int stbiw__zlib_bitrev(int code, int codebits)
{
   int res=0;
   while (codebits--) {
      res = (res << 1) | (code & 1);
      code >>= 1;
   }
   return res;
}

static unsigned int stbiw__zlib_countm(unsigned char *a, unsigned char *b, int limit)
{
   int i;
   for (i=0; i < limit && i < 258; ++i)
      if (a[i] != b[i]) break;
   return i;
}

static unsigned int stbiw__zhash(unsigned char *data)
{
   stbiw_uint32 hash = data[0] + (data[1] << 8) + (data[2] << 16);
   hash ^= hash << 3;
   hash += hash >> 5;
   hash ^= hash << 4;
   hash += hash >> 17;
   hash ^= hash << 25;
   hash += hash >> 6;
   return hash;
}

#define stbiw__zlib_flush() (out = stbiw__zlib_flushf(out, &bitbuf, &bitcount))
#define stbiw__zlib_add(code,codebits) \
      (bitbuf |= (code) << bitcount, bitcount += (codebits), stbiw__zlib_flush())
#define stbiw__zlib_huffa(b,c)  stbiw__zlib_add(stbiw__zlib_bitrev(b,c),c)
// default huffman tables
#define stbiw__zlib_huff1(n)  stbiw__zlib_huffa(0x30 + (n), 8)
#define stbiw__zlib_huff2(n)  stbiw__zlib_huffa(0x190 + (n)-144, 9)
#define stbiw__zlib_huff3(n)  stbiw__zlib_huffa(0 + (n)-256,7)
#define stbiw__zlib_huff4(n)  stbiw__zlib_huffa(0xc0 + (n)-280,8)
#define stbiw__zlib_huff(n)  ((n) <= 143 ? stbiw__zlib_huff1(n) : (n) <= 255 ? stbiw__zlib_huff2(n) : (n) <= 279 ? stbiw__zlib_huff3(n) : stbiw__zlib_huff4(n))
#define stbiw__zlib_huffb(n) ((n) <= 143 ? stbiw__zlib_huff1(n) : stbiw__zlib_huff2(n))

#define stbiw__ZHASH   16384

unsigned char * stbi_zlib_compress(unsigned char *data, int data_len, int *out_len, int quality)
{
   static unsigned short lengthc[] = { 3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258, 259 };
   static unsigned char  lengtheb[]= { 0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0 };
   static unsigned short distc[]   = { 1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577, 32768 };
   static unsigned char  disteb[]  = { 0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13 };
   unsigned int bitbuf=0;
   int i,j, bitcount=0;
   unsigned char *out = NULL;
   unsigned char **hash_table[stbiw__ZHASH]; // 64KB on the stack!
   if (quality < 5) quality = 5;

   stbiw__sbpush(out, 0x78);   // DEFLATE 32K window
   stbiw__sbpush(out, 0x5e);   // FLEVEL = 1
   stbiw__zlib_add(1,1);  // BFINAL = 1
   stbiw__zlib_add(1,2);  // BTYPE = 1 -- fixed huffman

   for (i=0; i < stbiw__ZHASH; ++i)
      hash_table[i] = NULL;

   i=0;
   while (i < data_len-3) {
      // hash next 3 bytes of data to be compressed
      int h = stbiw__zhash(data+i)&(stbiw__ZHASH-1), best=3;
      unsigned char *bestloc = 0;
      unsigned char **hlist = hash_table[h];
      int n = stbiw__sbcount(hlist);
      for (j=0; j < n; ++j) {
         if (hlist[j]-data > i-32768) { // if entry lies within window
            int d = stbiw__zlib_countm(hlist[j], data+i, data_len-i);
            if (d >= best) best=d,bestloc=hlist[j];
         }
      }
      // when hash table entry is too long, delete half the entries
      if (hash_table[h] && stbiw__sbn(hash_table[h]) == 2*quality) {
         STBIW_MEMMOVE(hash_table[h], hash_table[h]+quality, sizeof(hash_table[h][0])*quality);
         stbiw__sbn(hash_table[h]) = quality;
      }
      stbiw__sbpush(hash_table[h],data+i);

      if (bestloc) {
         // "lazy matching" - check match at *next* byte, and if it's better, do cur byte as literal
         h = stbiw__zhash(data+i+1)&(stbiw__ZHASH-1);
         hlist = hash_table[h];
         n = stbiw__sbcount(hlist);
         for (j=0; j < n; ++j) {
            if (hlist[j]-data > i-32767) {
               int e = stbiw__zlib_countm(hlist[j], data+i+1, data_len-i-1);
               if (e > best) { // if next match is better, bail on current match
                  bestloc = NULL;
                  break;
               }
            }
         }
      }

      if (bestloc) {
         int d = (int) (data+i - bestloc); // distance back
         STBIW_ASSERT(d <= 32767 && best <= 258);
         for (j=0; best > lengthc[j+1]-1; ++j);
         stbiw__zlib_huff(j+257);
         if (lengtheb[j]) stbiw__zlib_add(best - lengthc[j], lengtheb[j]);
         for (j=0; d > distc[j+1]-1; ++j);
         stbiw__zlib_add(stbiw__zlib_bitrev(j,5),5);
         if (disteb[j]) stbiw__zlib_add(d - distc[j], disteb[j]);
         i += best;
      } else {
         stbiw__zlib_huffb(data[i]);
         ++i;
      }
   }
   // write out final bytes
   for (;i < data_len; ++i)
      stbiw__zlib_huffb(data[i]);
   stbiw__zlib_huff(256); // end of block
   // pad with 0 bits to byte boundary
   while (bitcount)
      stbiw__zlib_add(0,1);

   for (i=0; i < stbiw__ZHASH; ++i)
      (void) stbiw__sbfree(hash_table[i]);

   {
      // compute adler32 on input
      unsigned int i=0, s1=1, s2=0, blocklen = data_len % 5552;
      int j=0;
      while (j < data_len) {
         for (i=0; i < blocklen; ++i) s1 += data[j+i], s2 += s1;
         s1 %= 65521, s2 %= 65521;
         j += blocklen;
         blocklen = 5552;
      }
      stbiw__sbpush(out, (unsigned char) (s2 >> 8));
      stbiw__sbpush(out, (unsigned char) s2);
      stbiw__sbpush(out, (unsigned char) (s1 >> 8));
      stbiw__sbpush(out, (unsigned char) s1);
   }
   *out_len = stbiw__sbn(out);
   // make returned pointer freeable
   STBIW_MEMMOVE(stbiw__sbraw(out), out, *out_len);
   return (unsigned char *) stbiw__sbraw(out);
}

unsigned int stbiw__crc32(unsigned char *buffer, int len)
{
   static unsigned int crc_table[256];
   unsigned int crc = ~0u;
   int i,j;
   if (crc_table[1] == 0)
      for(i=0; i < 256; i++)
         for (crc_table[i]=i, j=0; j < 8; ++j)
            crc_table[i] = (crc_table[i] >> 1) ^ (crc_table[i] & 1 ? 0xedb88320 : 0);
   for (i=0; i < len; ++i)
      crc = (crc >> 8) ^ crc_table[buffer[i] ^ (crc & 0xff)];
   return ~crc;
}

#define stbiw__wpng4(o,a,b,c,d) ((o)[0]=(unsigned char)(a),(o)[1]=(unsigned char)(b),(o)[2]=(unsigned char)(c),(o)[3]=(unsigned char)(d),(o)+=4)
#define stbiw__wp32(data,v) stbiw__wpng4(data, (v)>>24,(v)>>16,(v)>>8,(v));
#define stbiw__wptag(data,s) stbiw__wpng4(data, s[0],s[1],s[2],s[3])

static void stbiw__wpcrc(unsigned char **data, int len)
{
   unsigned int crc = stbiw__crc32(*data - len - 4, len+4);
   stbiw__wp32(*data, crc);
}

static unsigned char stbiw__paeth(int a, int b, int c)
{
   int p = a + b - c, pa = abs(p-a), pb = abs(p-b), pc = abs(p-c);
   if (pa <= pb && pa <= pc) return (unsigned char) a;
   if (pb <= pc) return (unsigned char) b;
   return (unsigned char) c;
}

unsigned char *stbi_write_png_to_mem(const unsigned char *pixels, int stride_bytes, int x, int y, int n, int *out_len)
{
   int ctype[5] = { -1, 0, 4, 2, 6 };
   unsigned char sig[8] = { 137,80,78,71,13,10,26,10 };
   unsigned char *out,*o, *filt, *zlib;
   signed char *line_buffer;
   int i,j,k,p,zlen;

   if (stride_bytes == 0)
      stride_bytes = x * n;

   filt = (unsigned char *) STBIW_MALLOC((x*n+1) * y); if (!filt) return 0;
   line_buffer = (signed char *) STBIW_MALLOC(x * n); if (!line_buffer) { STBIW_FREE(filt); return 0; }
   for (j=0; j < y; ++j) {
      static int mapping[] = { 0,1,2,3,4 };
      static int firstmap[] = { 0,1,0,5,6 };
      int *mymap = j ? mapping : firstmap;
      int best = 0, bestval = 0x7fffffff;
      for (p=0; p < 2; ++p) {
         for (k= p?best:0; k < 5; ++k) {
            int type = mymap[k],est=0;
            const unsigned char *z = pixels + stride_bytes*j;
            for (i=0; i < n; ++i)
               switch (type) {
                  case 0: line_buffer[i] = z[i]; break;
                  case 1: line_buffer[i] = z[i]; break;
                  case 2: line_buffer[i] = z[i] - z[i-stride_bytes]; break;
                  case 3: line_buffer[i] = z[i] - (z[i-stride_bytes]>>1); break;
                  case 4: line_buffer[i] = (signed char) (z[i] - stbiw__paeth(0,z[i-stride_bytes],0)); break;
                  case 5: line_buffer[i] = z[i]; break;
                  case 6: line_buffer[i] = z[i]; break;
               }
            for (i=n; i < x*n; ++i) {
               switch (type) {
                  case 0: line_buffer[i] = z[i]; break;
                  case 1: line_buffer[i] = z[i] - z[i-n]; break;
                  case 2: line_buffer[i] = z[i] - z[i-stride_bytes]; break;
                  case 3: line_buffer[i] = z[i] - ((z[i-n] + z[i-stride_bytes])>>1); break;
                  case 4: line_buffer[i] = z[i] - stbiw__paeth(z[i-n], z[i-stride_bytes], z[i-stride_bytes-n]); break;
                  case 5: line_buffer[i] = z[i] - (z[i-n]>>1); break;
                  case 6: line_buffer[i] = z[i] - stbiw__paeth(z[i-n], 0,0); break;
               }
            }
            if (p) break;
            for (i=0; i < x*n; ++i)
               est += abs((signed char) line_buffer[i]);
            if (est < bestval) { bestval = est; best = k; }
         }
      }
      // when we get here, best contains the filter type, and line_buffer contains the data
      filt[j*(x*n+1)] = (unsigned char) best;
      STBIW_MEMMOVE(filt+j*(x*n+1)+1, line_buffer, x*n);
   }
   STBIW_FREE(line_buffer);
   zlib = stbi_zlib_compress(filt, y*( x*n+1), &zlen, 8); // increase 8 to get smaller but use more memory
   STBIW_FREE(filt);
   if (!zlib) return 0;

   // each tag requires 12 bytes of overhead
   out = (unsigned char *) STBIW_MALLOC(8 + 12+13 + 12+zlen + 12);
   if (!out) return 0;
   *out_len = 8 + 12+13 + 12+zlen + 12;

   o=out;
   STBIW_MEMMOVE(o,sig,8); o+= 8;
   stbiw__wp32(o, 13); // header length
   stbiw__wptag(o, "IHDR");
   stbiw__wp32(o, x);
   stbiw__wp32(o, y);
   *o++ = 8;
   *o++ = (unsigned char) ctype[n];
   *o++ = 0;
   *o++ = 0;
   *o++ = 0;
   stbiw__wpcrc(&o,13);

   stbiw__wp32(o, zlen);
   stbiw__wptag(o, "IDAT");
   STBIW_MEMMOVE(o, zlib, zlen);
   o += zlen;
   STBIW_FREE(zlib);
   stbiw__wpcrc(&o, zlen);

   stbiw__wp32(o,0);
   stbiw__wptag(o, "IEND");
   stbiw__wpcrc(&o,0);

   STBIW_ASSERT(o == out + *out_len);

   return out;
}

int stbi_write_png(char const *filename, int x, int y, int comp, const void *data, int stride_bytes)
{
   FILE *f;
   int len;
   unsigned char *png = stbi_write_png_to_mem((const unsigned char *) data, stride_bytes, x, y, comp, &len);
   if (!png) return 0;
   f = fopen(filename, "wb");
   if (!f) { STBIW_FREE(png); return 0; }
   fwrite(png, 1, len, f);
   fclose(f);
   STBIW_FREE(png);
   return 1;
}

/* Revision history
      0.98 (2015-04-08)
             added STBIW_MALLOC, STBIW_ASSERT etc
      0.97 (2015-01-18)
             fixed HDR asserts, rewrote HDR rle logic
      0.96 (2015-01-17)
             add HDR output
             fix monochrome BMP
      0.95 (2014-08-17)
		       add monochrome TGA output
      0.94 (2014-05-31)
             rename private functions to avoid conflicts with stb_image.h
      0.93 (2014-05-27)
             warning fixes
      0.92 (2010-08-01)
             casts to unsigned char to fix warnings
      0.91 (2010-07-17)
             first public release
      0.90   first internal release
*/

#endif // Model_h
