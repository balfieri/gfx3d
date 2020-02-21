// Copyright (c) 2017-2018 Robert A. Alfieri
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
// Model.h - 3D model loading with one .h file
//
// How to use it:
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
//     4) If you want Model to generate compressed textures for you, it will use ARM's astcenc to encode
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
//     5) If you want Model to generate a BVH tree for you, add Model::BVH_TREE::BINARY as the third argument
//        to the Model() constructor in (1) to get a binary BVH tree.  QUAD and OCT trees are not
//        currently supported.
//
//        Note: BVH building reorders the polygons, so any polygon offsets in an Object structure will be garbage.
//              Most apps don't use the objects[] so this is not normally a problem.
//
// How it works:
//
//     1) Allocate large virtual memory 1D arrays for materials, texels, positions, normals, vertexes, polygons.
//        These are allocated on a page boundary to make uncompressed writes faster.
//        These arrays are dynamically resized.
//     2) Read entire .obj file into memory (the o/s should effectively make this work like an mmap).
//     3) Parse .obj file using custom parser that goes character-by-character and does its own number conversions. 
//     4) Add elements to 1D arrays.  Load any .mtl file encounted in .obj file.  
//     5) Optionally generate mipmap textures.  
//     6) Optionally compress textures.
//     7) Optionally generate BVH tree.
//     8) Write to  uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//     9) Read from uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//        The O/S will do the equivalent of an mmap() for each array.
//
// To Do:
//
//     0) Add some example programs to illustrate and test features.
//     1) ASTC texture decompression needs a lot more testing against ARM's reference model.
//     2) Support the .fbx binary and text formats.
//     3) Perhaps support other formats such as UE4 and MAX.
//     4) Test against more complex models, such as Disney Moana and Disney Cloud open-source models.
//     5) Add support for curves, patches, hair, etc.
//     6) Add a simple ray-tracing renderer as a separate file for an example and test.
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
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#endif

class Model
{
public:
    typedef int32_t  _int;                  // by default, we use 32-bit signed integers
    typedef uint32_t uint;                  // by default, we use 32-bit indexes for geometry
    typedef uint64_t uint64;                // by default, we use 64-bit indexes for texels
    typedef float    real;

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

    Model( std::string top_file, 
           MIPMAP_FILTER        mipmap_filter=MIPMAP_FILTER::NONE, 
           TEXTURE_COMPRESSION  texture_compression=TEXTURE_COMPRESSION::NONE,
           BVH_TREE             bvh_tree=BVH_TREE::NONE, 
           bool                 resolve_models=true );
    Model( std::string model_file, bool is_compressed );
    ~Model(); 

    bool write( std::string file_path, bool is_compressed ); 
    bool replace_materials( std::string mtl_file_path );

    // start of main structures
    static const uint VERSION = 0xB0BA1f09; // current version 

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    class real4
    {
    public:
        real c[4];
        
        real4( void )                               { c[0] = 0;  c[1] = 0;  c[2] = 0;  c[3] = 0;  }
        real4( real c0, real c1, real c2, real c3 ) { c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3; }

        real   dot( const real4 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const ;
        real4& normalize( void );
        real4  normalized( void ) const;
        real4  operator + ( const real4& v ) const;
        real4  operator - ( const real4& v ) const;
        real4  operator * ( const real4& v ) const;
        real4  operator * ( real s ) const;
        real4  operator / ( const real4& v ) const;
        real4  operator / ( real s ) const;
        bool   operator == ( const real4 &v2 ) const;
        real4& operator += ( const real4 &v2 );
        real4& operator -= ( const real4 &v2 );
        real4& operator *= ( const real4 &v2 );
        real4& operator *= ( const real s );
        real4& operator /= ( const real4 &v2 );
        real4& operator /= ( const real s );
    };

    class real3
    {
    public:
        real c[3];
        
        real3( void )                      { c[0] = 0;  c[1] = 0;  c[2] = 0;  }
        real3( real c0, real c1, real c2 ) { c[0] = c0; c[1] = c1; c[2] = c2; }

        real   dot( const real3 &v2 ) const;
        real3  cross( const real3 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const ;
        real3& normalize( void );
        real3  normalized( void ) const;
        real3  operator + ( const real3& v ) const;
        real3  operator - ( const real3& v ) const;
        real3  operator * ( const real3& v ) const;
        real3  operator * ( real s ) const;
        real3  operator / ( const real3& v ) const;
        real3  operator / ( real s ) const;
        bool   operator == ( const real3 &v2 ) const;
        real3& operator += ( const real3 &v2 );
        real3& operator -= ( const real3 &v2 );
        real3& operator *= ( const real3 &v2 );
        real3& operator *= ( const real s );
        real3& operator /= ( const real3 &v2 );
        real3& operator /= ( const real s );
    };

    class real2
    {
    public:
        real c[2];

        real2( void )             { c[0] = 0;  c[1] = 0;  }
        real2( real c0, real c1 ) { c[0] = c0; c[1] = c1; }

        real   dot( const real2 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const ;
        real2& normalize( void );
        real2  normalized( void ) const;
        real2  operator + ( const real2& v ) const;
        real2  operator - ( const real2& v ) const;
        real2  operator * ( const real2& v ) const;
        real2  operator * ( real s ) const;
        real2  operator / ( const real2& v ) const;
        real2  operator / ( real s ) const;
        real2& operator += ( const real2 &v2 );
        real2& operator -= ( const real2 &v2 );
        real2& operator *= ( const real2 &v2 );
        real2& operator *= ( const real s );
        real2& operator /= ( const real2 &v2 );
        real2& operator /= ( const real s );
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
        uint64      texel_cnt;              // in texels array  (last in file)

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
        real        opacity_scale;          // not sure what this is for
        real        shadow_caster_count;    // not sure what this is for
        real        tone_white;             // tone mapping white parameter
        real        tone_key;               // tone mapping key parameter

        uint64      camera_cnt;             // in cameras array
        uint64      initial_camera_i;       // initial active camera index (default: 0)
        uint64      frame_cnt;              // in frames array 
        uint64      animation_cnt;          // in animations array
        real        animation_speed;        // divide frame time by this to get real time (default: 1.0)

        // move these up later
        TEXTURE_COMPRESSION texture_compression; // NONE or ASTC used (individual textures are also marked)
    };

    // TODO: move these into HDR later
        real        tone_avg_luminance;     // tone mapping average luminance override (if supported by application)
        uint        direct_V;               // number of virtual point lights (VPLs) for direct lighting reservoir sampling (default: 0 == choose)
        uint        direct_M;               // number of M candidate samples for direct lighting reservoir sampling (default: 0 == choose)
        uint        direct_N;               // number of N used samples for direct lighting reservoir sampling (default: 0 == choose)

        uint64      graph_node_cnt;         // in graph_nodes array
        uint64      graph_node_cnt_max;     // in graph_nodes array
        uint64      graph_root_i;           // index of root Graph_Node in graph_nodes array

    class Object
    {
    public:
        uint        name_i;                 // index of object name in strings array
        uint        poly_cnt;               // number of polygons in this object
        uint        poly_i;                 // index of first polygon in polygons array
    };

    class HitInfo
    {
    public:
        uint            poly_i;
        real            t;  
        real            u;
        real            v;
        real            frac_uv_cov;            // UV footprint of hit surface (used by mip_level calculation)
        real3           p;
        real3           normal;
        real3           shading_normal;
        real3           tangent;
        real3           tangent_normalized;
        real3           bitangent;
        const Model *   model;              // model of owning BVH
    };

    class AABB                              // axis aligned bounding box
    {
    public:
        real3           min;                // bounding box min
        real3           max;                // bounding box max

        AABB( void ) {}
        AABB( const real3& p );             // init with one point
        AABB( const real3& p0, const real3& p1, const real3& p2 );

        void pad( real p );
        void expand( const AABB& other );
        void expand( const real3& p );
        bool encloses( const AABB& other ) const;
        bool hit( const real3& origin, const real3& direction, const real3& direction_inv, real tmin, real tmax ) const; 
    };

    class Polygon
    {
    public:
        uint        mtl_i;                  // index into materials array
        uint        vtx_cnt;                // number of vertices
        uint        vtx_i;                  // index into vertexes array of first vertex
        real3       normal;                 // surface normal
        real        area;                   // surface area };
        
        bool bounding_box( const Model * model, AABB& box, real padding=0 ) const;
        bool hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                  real solid_angle, real t_min, real t_max, HitInfo& hit_info ) const;
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

        inline bool is_emissive( void ) { return Ke.c[0] != 0.0 || Ke.c[1] != 0.0 || Ke.c[2] != 0.0; }
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

    class Texture
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
      //TEXTURE_COMPRESSION compression;    // NONE or ASTC (for now, look in model structure)
        uint            width;              // level 0 width
        uint            height;             // level 0 height
        uint            nchan;              // number of channels (typically 3 or 4)
        uint64          texel_i;            // index into texels array of first texel (uncompressed) or ASTC header
        real            bump_multiplier;    // should be for bump maps only but also used with other textures
        real3           offset;             // uvw offset of texture map on surface (w is 0 for 2D textures)
        real3           scale;              // uvw scale  of texture map on surface (w is 0 for 2D textures)

        real4           texel_read(              const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr=nullptr, uint64 * byte_cnt=nullptr,
                                                 bool do_srgb_to_linear=false ) const;
        static void     astc_blk_dims_get( uint width, uint height, uint& blk_dim_x, uint& blk_dim_y );

        // self-contained functions for implementing mip-based addressing 
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

    private:
        real4           texel_read_uncompressed( const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const;
        real4           texel_read_astc(         const Model * model, uint mip_level, uint64 ui, uint64 vi, 
                                                 uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const;

        real4           astc_decode( const unsigned char * bdata, const ASTC_Header * astc_hdr, uint s, uint t, uint r, bool do_srgb_to_linear ) const;
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
                  real solid_angle, real t_min, real t_max, HitInfo& hit_info ) const;
    };

    class Matrix                            // 4x4 transformation matrix used by Instance
    {
    public:
        real            m[4][4];

        Matrix( void )          { identity(); }

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
        void   transform( const real3& v, real3& r, bool div_by_w=true ) const; // r = *this * v (and divide by w by default)
        void   transform( const Matrix& M2, Matrix& M3 ) const; // M3 = *this * M2
        void   transpose( Matrix& mt ) const;           // return the transpose this matrix 
        void   invert( Matrix& minv ) const;            // return the inversion this matrix
        void   invert_affine( Matrix& minv ) const;     // same but faster because assumes an affine transform
        void   adjoint( Matrix& M ) const;              // used by invert()
        double determinant( void ) const;               // returns the determinant (as double for high-precision) 
        void   cofactor( Matrix& C ) const;             // used by adjoint() and determinant() 
        double subdeterminant( uint exclude_row, uint exclude_col ) const; // used by cofactor() 
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
                  real solid_angle, real t_min, real t_max, HitInfo& hit_info );
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
    };

    class Camera
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        real3           lookfrom;           // camera location
        real3           lookat;             // point camera is aimed at
        real3           vup;                // camera up direction
        real            aperture;           // lens aperture
        real            focus_dist;         // aka focal_length
        real            near;               // near clip plane distance (typical: 0.1)
        real            far;                // far clip plane distance (typical: 10000)
        real            vfov;               // vertical field of view
        real            aspect_ratio;       // aspect ratio
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

    uint string_alloc( std::string s );     // returns index in strings[] array

    // srgb<->linear conversion utilities
    static real srgb_to_linear_gamma( real srgb );
    static real srgb8_to_linear_gamma( uint8_t srgb );   // more efficient, uses LUT
    static real linear_to_srgb_gamma( real linear );
    static real linear8_to_srgb_gamma( uint8_t linear ); // more efficient, uses LUT

    // system utilities
    void cmd( std::string c, std::string error="command failed");  // calls std::system and aborts if not success

    // file utilities
    void dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); 
    bool file_exists( std::string file_name );
    bool file_read( std::string file_name, char *& start, char *& end );                        // sucks in entire file
    void file_write( std::string file_name, const unsigned char * data, uint64_t byte_cnt );

    // parsing utilities for files sucked into memory
    uint line_num;
    bool skip_whitespace_to_eol( char *& xxx, char *& xxx_end );  // on this line only
    bool skip_whitespace( char *& xxx, char *& xxx_end );
    bool skip_to_eol( char *& xxx, char *& xxx_end );
    bool eol( char *& xxx, char *& xxx_end );
    bool expect_char( char ch, char *& xxx, char* xxx_end, bool skip_whitespace_first=false );
    bool expect_cmd( const char * s, char *& xxx, char *& xxx_end );
    bool parse_string( std::string& s, char *& xxx, char *& xxx_end );
    bool parse_string_i( uint& s, char *& xxx, char *& xxx_end );
    bool parse_name( char *& name, char *& xxx, char *& xxx_end );
    bool parse_id( std::string& id, char *& xxx, char *& xxx_end );
    bool parse_id_i( uint& id_i, char *& xxx, char *& xxx_end );
    bool parse_option_name( std::string& option_name, char *& xxx, char *& xxx_end );
    bool parse_real3( real3& r3, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real2( real2& r2, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real( real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first=false );
    bool parse_int( _int& i, char *& xxx, char *& xxx_end );
    bool parse_uint( uint& u, char *& xxx, char *& xxx_end );
    bool parse_bool( bool& b, char *& xxx, char *& xxx_end );
    std::string surrounding_lines( char *& xxx, char *& xxx_end );

    // structs
    char *              mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
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
        CMD_MAP_BUMP
    } mtl_cmd_t;
            
    char * fsc_start;
    char * fsc_end;
    char * fsc;
    char * obj_start;
    char * obj_end;
    char * obj;
    char * mtl_start;
    char * mtl_end;
    char * mtl;
    char * gph_start;
    char * gph_end;
    char * gph;

    bool load_fsc( std::string fsc_file, std::string dir_name );        // .fscene 
    bool load_obj( std::string obj_file, std::string dir_name );        // .obj
    bool load_mtl( std::string mtl_file, std::string dir_name, bool replacing=false );
    bool load_tex( const char * tex_name, std::string dir_name, Texture *& texture );
    bool load_gph( std::string gph_file, std::string dir_name, std::string base_name, std::string ext_name );

    bool parse_obj_cmd( obj_cmd_t& cmd );
    bool parse_mtl_cmd( mtl_cmd_t& cmd, char *& mtl, char *& mtl_end );
    bool parse_gph_real( uint& r_i, char *& gph, char *& gph_end );
    bool parse_gph_real3( uint& r3_i, char *& gph, char *& gph_end );
    bool parse_gph_str( uint& s_i, char *& gph, char *& gph_end );
    bool parse_gph_id( uint& id_i, char *& gph, char *& gph_end );
    bool parse_gph_op_lookahead( uint& s_i, char *& gph, char *& gph_end );
    bool parse_gph_op( uint& op_i, char *& gph, char *& gph_end );
    bool parse_gph_term( uint& term_i, char *& gph, char *& gph_end );
    bool parse_gph_expr( uint& expr_i, char *& gph, char *& gph_end, uint curr_prec=0 );

    // BVH builder
    void bvh_build( BVH_TREE bvh_tree );
    uint bvh_qsplit( bool for_polys, uint poly_i, uint n, real pivot, uint axis );
    uint bvh_node( bool for_polys, uint i, uint n, uint axis );

    // allocates an array of T on a page boundary
    template<typename T>
    T * aligned_alloc( uint64 cnt );

    // reallocate array if we are about to exceed its current size
    template<typename T>
    inline void perhaps_realloc( T *& array, const uint64& hdr_cnt, uint64& max_cnt, uint64 add_cnt );

    bool write_uncompressed( std::string file_path );
    bool read_uncompressed( std::string file_path );
};

#define dprint( msg )
//#define dprint( msg ) std::cout << (msg) << "\n"
bool Model::debug = false;
#define mdout if (Model::debug) std::cout 

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n"; assert( false ); return false; }
#define fsc_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n"; assert( false ); goto error;   }
#define obj_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n"; assert( false ); goto error;   }
#define gph_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << "ERROR: " << msg << "\n"; assert( false ); goto error;   }
#define die_assert( bool, msg ) if ( !(bool) ) {                               std::cout << "ERROR: " << msg << "\n"; assert( false );               }
#define die( msg )                             {                               std::cout << "ERROR: " << msg << "\n"; assert( false );               }

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

inline std::string str( uint u )
{
    return std::to_string( u );
}

inline std::string str( Model::real r ) 
{
    return std::to_string( r );
}

inline std::string str( Model::real3 v ) 
{
    return "[" + str(v.c[0]) + "," + str(v.c[1]) + "," + str(v.c[2]) + "]";
}

inline std::ostream& operator << ( std::ostream& os, const Model::real3& v ) 
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
    os << box.min << ".." << box.max;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::Polygon& poly ) 
{
    os << "mtl_i=" << poly.mtl_i << " vtx_cnt=" << poly.vtx_cnt << " vtx_i=" << poly.vtx_i <<
          " normal=" << poly.normal << " area=" << poly.area;
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
          " map_Bump_i=" << mat.map_Bump_i;
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
    uint w = (hdr.xsize[0] << 0) | (hdr.xsize[1] << 8) | (hdr.xsize[2] << 16);
    uint h = (hdr.ysize[0] << 0) | (hdr.ysize[1] << 8) | (hdr.ysize[2] << 16);
    uint d = (hdr.zsize[0] << 0) | (hdr.zsize[1] << 8) | (hdr.zsize[2] << 16);
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

Model::Model( std::string                top_file,      
              Model::MIPMAP_FILTER       mipmap_filter, 
              Model::TEXTURE_COMPRESSION texture_compression,
              Model::BVH_TREE bvh_tree,  bool resolve_models )
{
    is_good = false;
    error_msg = "<unknown error>";
    mapped_region = nullptr;
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
         graph_root_i = uint(-1);
    hdr->lighting_scale = 1.0;
    hdr->ambient_intensity = real3( 0.1, 0.1, 0.1 );
    hdr->sky_box_tex_i = uint(-1);
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
         graph_node_cnt_max= 1;
         graph_node_cnt    = 0;
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
    graph_nodes       = aligned_alloc<Graph_Node>( graph_node_cnt_max );
    textures          = aligned_alloc<Texture>(  max->tex_cnt );
    texels            = aligned_alloc<unsigned char>( max->texel_cnt );
    bvh_nodes         = aligned_alloc<BVH_Node>( max->bvh_node_cnt );
    matrixes          = aligned_alloc<Matrix>(   max->matrix_cnt );
    instances         = aligned_alloc<Instance>( max->inst_cnt );
    lights            = aligned_alloc<Light>(    max->light_cnt );
    cameras           = aligned_alloc<Camera>(   max->camera_cnt );
    frames            = aligned_alloc<Frame>(    max->frame_cnt );
    animations        = aligned_alloc<Animation>(max->animation_cnt );

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
    } else {
        error_msg = "unknown top file ext_name: " + ext_name;
        return;
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
                                std::cout << "Ignoring poly_i=" << poly_i << " in " << file_name << "\n";
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
                M->transform( t_box->min, instance->box.min );
                M->transform( t_box->max, instance->box.max );
            }
        }
    }

    //------------------------------------------------------------
    // Optionally build BVH around Instances and/or Polygons
    //------------------------------------------------------------
    if ( bvh_tree != BVH_TREE::NONE ) bvh_build( bvh_tree );

    //------------------------------------------------------------
    // Add up byte count.
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
                    uint64( graph_node_cnt            ) * sizeof( graph_nodes[0] ) +
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
    graph_node_cnt = 0;         // for now
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
        _read( graph_nodes,         Graph_Node,    graph_node_cnt );
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
    // Recreate name maps.
    //------------------------------------------------------------
    for ( uint i = 0; i < hdr->obj_cnt; i++ )
    {
        std::string name = &strings[objects[i].name_i];
        name_to_obj_i[name] = i;
        mdout << "obj_i=" << i << " name_i=" << objects[i].name_i << " name=" << name << "\n";
    }

    for ( uint i = 0; i < hdr->mtl_cnt; i++ )
    {
        std::string name = &strings[materials[i].name_i];
        name_to_mtl_i[name] = i;
        mdout << "mtl_i=" << i << " name_i=" << materials[i].name_i << " name=" << name << "\n";
    }

    for ( uint i = 0; i < hdr->tex_cnt; i++ )
    {
        std::string name = &strings[textures[i].name_i];
        name_to_tex_i[name] = i;
        mdout << "tex_i=" << i << " name_i=" << textures[i].name_i << " name=" << name << "\n";
    }

    is_good = true;
}

Model::~Model()
{
    if ( mapped_region != nullptr ) {
        delete mapped_region;
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
    _write( graph_nodes,   	graph_node_cnt  	* sizeof(graph_nodes[0]) );
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
    _uwrite( graph_nodes,   	graph_node_cnt  	* sizeof(graph_nodes[0]) );
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
    if ( !file_read( model_path, start, end ) ) return false;
    mapped_region = start;

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
    _uread( graph_nodes,         Graph_Node,    graph_node_cnt );
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
    tone_avg_luminance = 0.0;
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
                
                } else if ( strcmp( field, "env_map_intensity_scale" ) == 0 ) {
                    if ( !parse_real( hdr->env_map_intensity_scale, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "opacity_scale" ) == 0 ) {
                    if ( !parse_real( hdr->opacity_scale, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "shadow_caster_count" ) == 0 ) {
                    if ( !parse_real( hdr->shadow_caster_count, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "background_color" ) == 0 ) {
                    if ( !parse_real3( hdr->background_color, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "direct_V" ) == 0 ) {
                    if ( !parse_uint( direct_V, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "direct_M" ) == 0 ) {
                    if ( !parse_uint( direct_M, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "direct_N" ) == 0 ) {
                    if ( !parse_uint( direct_N, fsc, fsc_end ) ) goto error;

                } else if ( strcmp( field, "tone_white" ) == 0 ) {
                    if ( !parse_real( hdr->tone_white, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "tone_key" ) == 0 ) {
                    if ( !parse_real( hdr->tone_key, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "tone_avg_luminance" ) == 0 ) {
                    if ( !parse_real( tone_avg_luminance, fsc, fsc_end, true ) ) goto error;

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
                camera->aspect_ratio = 1.0;

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
                        if ( !parse_real( camera->aspect_ratio, fsc, fsc_end, true ) ) goto error;

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
                            camera->aspect_ratio = -1.0;

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
                                    if ( !parse_real( camera->aspect_ratio, fsc, fsc_end, true ) ) goto error;

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
    char *      mtllib = nullptr;
    char *      obj_name = nullptr;
    char *      name;
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

    char * mtl_name = nullptr;
    char * tex_name = nullptr;

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

bool Model::load_tex( const char * tex_name, std::string dir_name, Model::Texture *& texture )
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

    if ( hdr->texture_compression == TEXTURE_COMPRESSION::ASTC ) {
        //---------------------------------------------
        // Align to 16B boundary to match ASTC alignment
        // Should not need to reallocate due to >16B granularity of allocations.
        //---------------------------------------------
        hdr->texel_cnt += ((hdr->texel_cnt % 16) == 0) ? 0 : (16 - (hdr->texel_cnt % 16)); 
    }
    texture->texel_i = hdr->texel_cnt;  // always

    unsigned char * data = nullptr;
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
                int    w;
                int    h;
                int    nchan;
                data = stbi_load( file_name, &w, &h, &nchan, 0 );
                rtn_assert( data != nullptr, "unable to read in texture file " + std::string( file_name ) );
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
                                unsigned char * frgb = data + fii*byte_width + fjj*texture->nchan;
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
    if ( data != nullptr ) delete data;
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
    graph_root_i = block_i;

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

error:
    error_msg += " (at line " + std::to_string( line_num ) + " of " + gph_file + ")";
    die_assert( false, "aborting..." );
    return false;
}

uint Model::gph_node_alloc( Model::GRAPH_NODE_KIND kind, uint parent_i )
{
    //mdout << "new node: kind=" << kind << "\n";
    perhaps_realloc( graph_nodes, graph_node_cnt, graph_node_cnt_max, 1 );
    uint node_i = graph_node_cnt++;
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

bool Model::file_read( std::string file_path, char *& start, char *& end )
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

    // this large read should behave like an mmap() inside the o/s kernel and be as fast
    start = aligned_alloc<char>( size );
    if ( start == nullptr ) {
        close( fd );
        rtn_assert( 0, "could not read file " + std::string(fname) + " - malloc() error: " + strerror( errno ) );
    }
    end = start + size;

    char * addr = start;
    while( size != 0 ) 
    {
        size_t _this_size = 1024*1024*1024;
        if ( size < _this_size ) _this_size = size;
        if ( ::read( fd, addr, _this_size ) <= 0 ) {
            close( fd );
            rtn_assert( 0, "could not read() file " + std::string(fname) + " - read error: " + strerror( errno ) );
        }
        size -= _this_size;
        addr += _this_size;
    }
    close( fd );
    return true;
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

void Model::cmd( std::string c, std::string error )
{
    std::cout << c << "\n";
    if ( std::system( c.c_str() ) != 0 ) die_assert( false, "ERROR: " + error + ": " + c );
}

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

inline bool Model::skip_whitespace( char *& xxx, char *& xxx_end )
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
                for( char * ccc = xxx+1; ccc != xxx_end && *ccc != '\n' && *ccc != '\r'; ccc++ )
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

inline bool Model::skip_whitespace_to_eol( char *& xxx, char *& xxx_end )
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

inline bool Model::skip_to_eol( char *& xxx, char *& xxx_end )
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

inline bool Model::eol( char *& xxx, char *& xxx_end )
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

inline bool Model::expect_char( char ch, char *& xxx, char* xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );
    rtn_assert( xxx != xxx_end, "premature end of file" );
    rtn_assert( *xxx == ch, "expected character '" + std::string(1, ch) + "' got '" + std::string( 1, *xxx ) + "' " + surrounding_lines( xxx, xxx_end ) );
    xxx++;
    return true;
}

inline bool Model::expect_cmd( const char * s, char *& xxx, char *& xxx_end )
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

inline bool Model::parse_string( std::string& s, char *& xxx, char *& xxx_end )
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

inline uint Model::string_alloc( std::string s )
{
    uint s_len = s.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    uint s_i = hdr->char_cnt;
    char * to_s = &strings[hdr->char_cnt];
    hdr->char_cnt += s_len + 1;
    memcpy( to_s, s.c_str(), s_len+1 );
    return s_i;
}

inline bool Model::parse_string_i( uint& s_i, char *& xxx, char *& xxx_end )
{
    std::string s;
    if ( !parse_string( s, xxx, xxx_end ) ) return false;

    s_i = string_alloc( s );
    return true;
}

inline bool Model::parse_name( char *& name, char *& xxx, char *& xxx_end )
{
    bool vld = false;
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, 1024 );
    name = &strings[hdr->char_cnt];

    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    uint len = 0;
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( ch == '\n' || ch == '\r' ) break;

        rtn_assert( len < 1024, "string is larger than 1024 characters" );
        name[len++] = ch;
        vld = true;
        xxx++;
    }

    if ( vld ) {
        name[len] = '\0';
        hdr->char_cnt += len+1;
        char * ptr;
        for( ptr = &name[len-1]; ptr != name; ptr-- )
        {
            // skip trailing spaces
            if ( *ptr != ' ' && *ptr != '\t' ) break;
            *ptr = '\0';
        }

        return *ptr != '\0';
    }

    rtn_assert( 0, "could not parse name: " + surrounding_lines( xxx, xxx_end ) );
}

inline bool Model::parse_id( std::string& id, char *& xxx, char *& xxx_end )
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

inline bool Model::parse_id_i( uint& id_i, char *& xxx, char *& xxx_end )
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

inline bool Model::parse_option_name( std::string& option_name, char *& xxx, char *& xxx_end )
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

inline bool Model::parse_mtl_cmd( mtl_cmd_t& cmd, char *& mtl, char *& mtl_end )
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

        default:
            rtn_assert( 0, "bad .mtl command character '" + std::string( 1, ch ) + "': " + surrounding_lines( mtl, mtl_end ) );
    }
}

inline bool Model::parse_gph_real( uint& r_i, char *& gph, char *& gph_end )
{
    real r;
    if ( !parse_real( r, gph, gph_end, true ) ) return false;

    r_i = gph_node_alloc( GRAPH_NODE_KIND::REAL );
    mdout << "    " << r << "\n";
    graph_nodes[r_i].u.r = r;

    return true;
}

inline bool Model::parse_gph_real3( uint& r3_i, char *& gph, char *& gph_end )
{
    real3 r3;
    if ( !parse_real3( r3, gph, gph_end, true ) ) return false;

    r3_i = gph_node_alloc( GRAPH_NODE_KIND::REAL3 );
    mdout << "    " << r3 << "\n";
    graph_nodes[r3_i].u.r3 = r3;

    return true;
}

inline bool Model::parse_gph_str( uint& s_i, char *& gph, char *& gph_end )
{
    uint str_i;
    if ( !parse_string_i( str_i, gph, gph_end ) ) return false;

    s_i = gph_node_alloc( GRAPH_NODE_KIND::STR );
    mdout << "    " << &strings[str_i] << "\n";
    graph_nodes[s_i].u.s_i = str_i;

    return true;
}

inline bool Model::parse_gph_id( uint& id_i, char *& gph, char *& gph_end )
{
    uint s_i;
    if ( !parse_id_i( s_i, gph, gph_end ) ) return false;

    id_i = gph_node_alloc( GRAPH_NODE_KIND::ID );
    mdout << "    " << &strings[s_i] << "\n";
    graph_nodes[id_i].u.s_i = s_i;

    return true;
}

inline bool Model::parse_gph_op_lookahead( uint& s_i, char *& gph, char *& gph_end )
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

inline bool Model::parse_gph_op( uint& op_i, char *& gph, char *& gph_end )
{
    uint s_i;
    if ( !parse_gph_op_lookahead( s_i, gph, gph_end ) ) return false;
    gph += strlen( &strings[s_i] );

    op_i = gph_node_alloc( GRAPH_NODE_KIND::OP );
    graph_nodes[op_i].u.s_i = s_i;
    mdout << "    " << std::string(&strings[s_i]) << "\n";

    return true;
}

bool Model::parse_gph_term( uint& term_i, char *& gph, char *& gph_end )
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

inline bool Model::parse_gph_expr( uint& expr_i, char *& gph, char *& gph_end, uint curr_prec )
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

inline bool Model::parse_real3( Model::real3& r3, char *& xxx, char *& xxx_end, bool has_brackets )
{
    return (!has_brackets || expect_char( '[', xxx, xxx_end, true )) &&
           parse_real( r3.c[0], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[1], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[2], xxx, xxx_end, has_brackets ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Model::parse_real2( Model::real2& r2, char *& xxx, char *& xxx_end, bool has_brackets )
{
    return (!has_brackets || expect_char( '[', xxx, xxx_end, true )) &&
           parse_real( r2.c[0], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r2.c[1], xxx, xxx_end, has_brackets ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Model::parse_real( Model::real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first )
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
            //r = std::nan( "1" );
            r = 0.0;                    // make them zeros
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

    r = std::atof( s.c_str() );
    dprint( "real=" + std::to_string( r ) );
    return true;
}

inline bool Model::parse_int( _int& i, char *& xxx, char *& xxx_end )
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

inline bool Model::parse_uint( uint& u, char *& xxx, char *& xxx_end )
{
    _int i;
    if ( !parse_int( i, xxx, xxx_end ) ) return false;
    rtn_assert( i >= 0, "parse_uint encountered negative integer" );
    u = i;
    return true;
}

inline bool Model::parse_bool( bool& b, char *& xxx, char *& xxx_end )
{
    std::string id;
    if ( !parse_id( id , xxx, xxx_end ) ) return false;
    b = id == std::string( "true" );
    return true;
}

std::string Model::surrounding_lines( char *& xxx, char *& xxx_end )
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

inline Model::AABB::AABB( const Model::real3& p )
{
    min = p;
    max = p;
}  

inline Model::AABB::AABB( const Model::real3& p0, const Model::real3& p1, const Model::real3& p2 ) 
{
    min = p0;
    max = p0;
    expand( p1 );
    expand( p2 );
}  

void Model::bvh_build( Model::BVH_TREE bvh_tree )
{
    (void)bvh_tree;
    if ( hdr->poly_cnt != 0 ) {
        die_assert( hdr->inst_cnt == 0, "inst_cnt should be 0" );
        hdr->bvh_root_i = bvh_node( true, 0, hdr->poly_cnt, 1 );

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
    } else {
        die_assert( hdr->inst_cnt != 0 && hdr->poly_cnt == 0, "poly_cnt should be 0" );
        hdr->bvh_root_i = bvh_node( false, 0, hdr->inst_cnt, 1 );
    }
}

inline Model::uint Model::bvh_qsplit( bool for_polys, Model::uint first, Model::uint n, Model::real pivot, Model::uint axis )
{
   uint m = first;

   for( uint i = first; i < (first+n); i++ )
   {
       AABB box;
       polygons[i].bounding_box( this, box );
       real centroid = (box.min.c[axis] + box.max.c[axis]) * 0.5;
       if ( centroid < pivot ) {
           if ( for_polys ) {
               Polygon temp = polygons[i];
               polygons[i]  = polygons[m];
               polygons[m]  = temp;
           } else {
               Instance temp = instances[i];
               instances[i]  = instances[m];
               instances[m]  = temp;
           }
           m++;
       }
    }

    die_assert( m >= first && m <= (first+n), "qsplit has gone mad" );
    if ( m == first || m == (first +n) ) m = first + n/2;
    return m;
}

Model::uint Model::bvh_node( bool for_polys, Model::uint first, Model::uint n, Model::uint axis ) 
{
    perhaps_realloc( bvh_nodes, hdr->bvh_node_cnt, max->bvh_node_cnt, 1 );
    uint bvh_i = hdr->bvh_node_cnt++;
    BVH_Node * node = &bvh_nodes[bvh_i];

    if ( for_polys ) {
        polygons[first].bounding_box( this, node->box );
    } else {
        instances[first].bounding_box( this, node->box );
    }
    for( uint i = 1; i < n; i++ )
    {
        AABB new_box;
        if ( for_polys ) {
            polygons[first+i].bounding_box( this, new_box );
        } else {
            instances[first+i].bounding_box( this, new_box );
        }
        node->box.expand( new_box );
    }

    if ( n == 1 || n == 2 ) {
        node->left_i = first;
        AABB new_box;
        if ( for_polys ) {
            node->left_kind  = Model::BVH_NODE_KIND::POLYGON;
            node->right_kind = Model::BVH_NODE_KIND::POLYGON;
            polygons[first+0].bounding_box( this, new_box );
        } else {
            node->left_kind  = Model::BVH_NODE_KIND::INSTANCE;
            node->right_kind = Model::BVH_NODE_KIND::INSTANCE;
            instances[first+0].bounding_box( this, new_box );
        }
        die_assert( node->box.encloses( new_box ), "box should enclose new_box" );
        if ( n == 2 ) {
            if ( for_polys ) {
                polygons[first+1].bounding_box( this, new_box );
            } else {
                instances[first+1].bounding_box( this, new_box );
            }
            die_assert( node->box.encloses( new_box ), "box should enclose new_box" );
            node->right_i = first + 1;
        } else {
            node->right_i = first;
        }

    } else {
        node->left_kind = Model::BVH_NODE_KIND::BVH_NODE;
        node->right_kind = Model::BVH_NODE_KIND::BVH_NODE;
        real pivot = (node->box.min.c[axis] + node->box.max.c[axis]) * 0.5;
        uint m = bvh_qsplit( for_polys, first, n, pivot, axis );
        uint nm = m - first;
        uint left_i  = bvh_node( for_polys, first, nm,   (axis + 1) % 3 );
        uint right_i = bvh_node( for_polys, m,     n-nm, (axis + 1) % 3 );
        node = &bvh_nodes[bvh_i];  // could change after previous calls
        node->left_i  = left_i;
        node->right_i = right_i;
        die_assert( node->box.encloses( bvh_nodes[left_i].box ) && node->box.encloses( bvh_nodes[right_i].box ), "box does not enclose left and right boxes" );
    }

    return bvh_i;
}

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

inline Model::real4& Model::real4::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real4 Model::real4::normalized( void ) const
{
    return *this / length();
}

inline Model::real4 Model::real4::operator + ( const Model::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    r.c[2] = c[2] + v2.c[2];
    r.c[3] = c[3] + v2.c[3];
    return r;
}

inline Model::real4 Model::real4::operator - ( const Model::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    r.c[2] = c[2] - v2.c[2];
    r.c[3] = c[3] - v2.c[3];
    return r;
}

inline Model::real4 Model::real4::operator * ( const Model::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    r.c[2] = c[2] * v2.c[2];
    r.c[3] = c[3] * v2.c[3];
    return r;
}

inline Model::real4 operator * ( Model::real s, const Model::real4& v ) 
{
    return Model::real4( s*v.c[0], s*v.c[1], s*v.c[2], s*v.c[3] );
}

inline Model::real4 Model::real4::operator * ( Model::real s ) const
{
    real4 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    r.c[2] = c[2] * s;
    r.c[3] = c[3] * s;
    return r;
}

inline Model::real4 Model::real4::operator / ( const Model::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    r.c[2] = c[2] / v2.c[2];
    r.c[3] = c[3] / v2.c[3];
    return r;
}

inline Model::real4 Model::real4::operator / ( Model::real s ) const
{
    real4 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    r.c[2] = c[2] / s;
    r.c[3] = c[3] / s;
    return r;
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

inline Model::real3& Model::real3::normalize( void )
{
    *this /= length();
    return *this;
}

inline Model::real3 Model::real3::normalized( void ) const
{
    return *this / length();
}

inline Model::real3 Model::real3::operator + ( const Model::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    r.c[2] = c[2] + v2.c[2];
    return r;
}

inline Model::real3 Model::real3::operator - ( const Model::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    r.c[2] = c[2] - v2.c[2];
    return r;
}

inline Model::real3 Model::real3::operator * ( const Model::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    r.c[2] = c[2] * v2.c[2];
    return r;
}

inline Model::real3 operator * ( Model::real s, const Model::real3& v ) 
{
    return Model::real3( s*v.c[0], s*v.c[1], s*v.c[2] );
}

inline Model::real3 Model::real3::operator * ( Model::real s ) const
{
    real3 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    r.c[2] = c[2] * s;
    return r;
}

inline Model::real3 Model::real3::operator / ( const Model::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    r.c[2] = c[2] / v2.c[2];
    return r;
}

inline Model::real3 Model::real3::operator / ( Model::real s ) const
{
    real3 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    r.c[2] = c[2] / s;
    return r;
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

inline Model::real2 Model::real2::operator + ( const Model::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    return r;
}

inline Model::real2 Model::real2::operator - ( const Model::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    return r;
}

inline Model::real2 Model::real2::operator * ( const Model::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    return r;
}

inline Model::real2 Model::real2::operator * ( Model::real s ) const
{
    real2 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    return r;
}

inline Model::real2 Model::real2::operator / ( const Model::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    return r;
}

inline Model::real2 Model::real2::operator / ( Model::real s ) const
{
    real2 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    return r;
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

void Model::Matrix::transform( const real4& v, real4& r ) const
{
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

inline Model::real4 Model::Texture::texel_read( const Model * model, uint mip_level, uint64 ui, uint64 vi, uint64 * vaddr, uint64 * byte_cnt, bool do_srgb_to_linear ) const
{
    TEXTURE_COMPRESSION compression = model->hdr->texture_compression;
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
    real4 rgba = astc_decode( bdata, astc_hdr, s, t, r, do_srgb_to_linear );

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
                                        real      frac_uv_covg, 
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
    const uint64_t sqrt_frac_uv_covg_fxd = int( std::sqrt( frac_uv_covg ) * real(1 << covg_frac_w) );

    mdout << "astc_blk_addr_get: ray_id=" << ray_id <<
                               " blockdim_x=" << blockdim_x << " blockdim_y=" << blockdim_y << " blockdim_z=" << blockdim_z <<
                               std::hex << " xsize=0x" << xsize << " ysize=0x" << ysize << " zsize=0x" << zsize << 
                               " sqrt_xsize_ysize_zsize=0x" << sqrt_xsize_ysize_zsize << 
                               std::dec << " u=" << u << " v=" << v << " w=" << w << 
                               std::hex << " u_fxd=0x" << u_fxd << " v_fxd=0x" << v_fxd << " w_fxd=0x" << w_fxd << 
                               " frac_uv_covg=" << frac_uv_covg << " sqrt_frac_uv_covg_fxd=0x" << sqrt_frac_uv_covg_fxd <<
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
    const uint64_t iwidth_of_footprint_lg2 = uint( std::log2( real(iwidth_of_footprint) ) );
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

        mdout << "mip" << i << ": ray_id=" << ray_id << 
                std::hex << " mip_blk0_offset=0x" << mip_blk0_offset << " mip_blk_cnt_x=0x" << mip_blk_cnt_x << " mip_blk_cnt_y=0x" << mip_blk_cnt_y << " mip_blk_cnt_z=0x" << mip_blk_cnt_z << 
                " mip_blk_cnt=0x" << mip_blk_cnt << 
                std::dec << "\n";

        if ( i == (size_w-1) || i >= iwidth_of_footprint_lg2 ) break;

        mip_blk0_offset += mip_blk_cnt;
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
    uint64_t vi = (v_fxd * mip_ysize) >> uv_frac_w;
    uint64_t wi = (w_fxd * mip_zsize) >> uv_frac_w;
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

inline Model::real4 Model::Texture::astc_decode( const unsigned char * bdata, const ASTC_Header * astc_hdr, uint s, uint t, uint r, bool do_srgb_to_linear ) const
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
            die_assert( false, "ASTC HDR color endpoint modes are not yet supported" );
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
    die_assert( vi < cnt, "astc_decode_integer vi >= cnt" );
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

inline void Model::AABB::pad( Model::real p ) 
{
    min -= real3( p, p, p );
    max += real3( p, p, p );
}

inline void Model::AABB::expand( const Model::AABB& other )
{
    for( uint i = 0; i < 3; i++ )
    {
        if ( other.min.c[i] < min.c[i] ) min.c[i] = other.min.c[i];
        if ( other.max.c[i] > max.c[i] ) max.c[i] = other.max.c[i];
    }
}

inline void Model::AABB::expand( const Model::real3& p ) 
{
    if ( p.c[0] < min.c[0] ) min.c[0] = p.c[0];
    if ( p.c[1] < min.c[1] ) min.c[1] = p.c[1];
    if ( p.c[2] < min.c[2] ) min.c[2] = p.c[2];
    if ( p.c[0] > max.c[0] ) max.c[0] = p.c[0];
    if ( p.c[1] > max.c[1] ) max.c[1] = p.c[1];
    if ( p.c[2] > max.c[2] ) max.c[2] = p.c[2];
}

inline bool Model::AABB::encloses( const AABB& other ) const
{
    return min.c[0] <= other.min.c[0] &&
           min.c[1] <= other.min.c[1] &&
           min.c[2] <= other.min.c[2] &&
           max.c[0] >= other.max.c[0] &&
           max.c[1] >= other.max.c[1] &&
           max.c[2] >= other.max.c[2];
}

inline bool Model::AABB::hit( const Model::real3& origin, const Model::real3& direction, const Model::real3& direction_inv, 
                              Model::real tmin, Model::real tmax ) const 
{
    mdout << "Model::AABB::hit: " << *this << " tmin=" << tmin << " tmax=" << tmax << "\n";
    (void)direction;
    for( uint a = 0; a < 3; a++ ) 
    {
        real dir_inv = direction_inv.c[a];
        real v0 = (min.c[a] - origin.c[a]) * dir_inv;
        real v1 = (max.c[a] - origin.c[a]) * dir_inv;
        tmin = std::fmax( tmin, std::fmin( v0, v1 ) );
        tmax = std::fmin( tmax, std::fmax( v0, v1 ) );
        mdout << "Model::AABB::hit:     " << a << ": min=" << min.c[a] << " max=" << max.c[a] << 
                                   " dir_inv=" << dir_inv << " origin=" << origin.c[a] << 
                                   " v0=" << v0 << " v1=" << v1 << " tmin=" << tmin << " tmax=" << tmax << "\n";
    }
    bool r = tmax >= std::fmax( tmin, real(0.0) );
    mdout << "Model::AABB::hit: return=" << r << "\n";
    return r;
}

inline bool Model::Polygon::hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv,
        real solid_angle, real t_min, real t_max, HitInfo& hit_info ) const
{
    (void)direction_inv;
    if ( vtx_cnt == 3 ) {
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
        mdout << "Model::Polygon::hit: poly_i=" << poly_i << " origin=" << origin << 
                                     " direction=" << direction << " direction_inv=" << direction_inv << " solid_angle=" << solid_angle <<
                                     " d=" << d << " t=" << t << " t_min=" << t_min << " t_max=" << t_max << "\n";
        if ( t > t_min && t < t_max ) {
            // compute barycentrics, see if it's in triangle
            const real3& p1 = positions[vertexes[1].v_i];
            const real3& p2 = positions[vertexes[2].v_i];
            const real3 p = origin + direction*t;
            // careful about order!
            const real3 p_m_p0  = p  - p0;
            const real3 p1_m_p0 = p1 - p0;
            const real3 p2_m_p0 = p2 - p0;
            real area1 = p1_m_p0.cross( p_m_p0 ).dot( normal );
            real area2 = p_m_p0.cross( p2_m_p0 ).dot( normal );
            real area_2x = area * 2.0;
            real beta    = area1/area_2x;
            real gamma   = area2/area_2x;
            mdout << "Model::Polygon::hit: poly_i=" << poly_i << " beta=" << beta << " gamma=" << gamma << "\n";
            if ( beta >= 0.0 && gamma >= 0.0 && (beta + gamma) <= 1.0 ) {
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
                hit_info.poly_i = poly_i;
                hit_info.t = t;

                 hit_info.normal = normal;
                if ( vertexes[0].vn_i != uint(-1) ) {
                    hit_info.shading_normal = n0*alpha + n1*gamma + n2*beta;
                    hit_info.shading_normal.normalize();
                } else {
                    hit_info.shading_normal = normal;
                }
                    hit_info.normal = normal;
                // epsilon to get it off the polygon
                // may cause trouble in two-sided... move to shader?
            //    hit_info.p = p + 0.01*hit_info.normal;;
                hit_info.p = p;

                real distance_squared = t*t / direction.length_sqr();
                real ray_footprint_area_on_triangle = solid_angle*distance_squared;
                real twice_uv_area_of_triangle = std::abs(0.5*(u0*v1 + u1*v2 + u2*v0 - u0*v2 - u1*v0 - u2*v1));
                hit_info.frac_uv_cov = ray_footprint_area_on_triangle * twice_uv_area_of_triangle / (2.0*area);


                hit_info.u = alpha*u0 + gamma*u1 + beta*u2 ;
                hit_info.v = alpha*v0 + gamma*v1 + beta*v2 ;

                real3 deltaPos1 = p1-p0;
                real3 deltaPos2 = p2-p0;

                real deltaU1 = u1-u0;
                real deltaU2 = u2-u0;
                real deltaV1 = v1-v0;
                real deltaV2 = v2-v0;

                real r = 1.0 / (deltaU1 * deltaV2 - deltaV1 * deltaU2);
                hit_info.tangent = (deltaPos1 * deltaV2   - deltaPos2 * deltaV1)*r;
                hit_info.tangent_normalized = hit_info.tangent.normalize();
                hit_info.bitangent = (deltaPos2 * deltaU1   - deltaPos1 * deltaU2)*r;

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
                    real u = hit_info.u;
                    real v = hit_info.v;
                    real sqrt_nx_ny = std::sqrt(nx*ny);
                    real width_of_footprint = std::sqrt(hit_info.frac_uv_cov) * sqrt_nx_ny;
                    real mip_level = std::log2( width_of_footprint );
                    int nchan = model->textures->nchan;
                    for (int imip_level = mip_level; imip_level > 0 && !(mx == 1 && my == 1); imip_level--)
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
                    /*
                    if (i >= mx) i -= mx;
                    if (j >= my) j -= my;
                    */
                    while (i >= mx) i -= mx;
                    while (j >= my) j -= my;

                    float opacity = float(mdata[nchan*i + nchan*mx*j+3]) / 255.0;

                    if (float(uniform()) > opacity) {
                        return false;
                    }
                }

                hit_info.model = model;

                mdout << "Model::Polygon::hit: poly_i=" << poly_i << " HIT t=" << hit_info.t << 
                         " p=" << hit_info.p << " normal=" << hit_info.normal << 
                         " frac_uv_cov=" << hit_info.frac_uv_cov << 
                         " u=" << hit_info.u << " v=" << hit_info.v << " mtl_i=" << mtl_i << "\n";
                return true;
            }
        }
    }
    mdout << "Model::Polygon::hit: NOT a hit, poly_i=" << (this - model->polygons) << "\n";
    return false;
}

bool Model::Instance::bounding_box( const Model * model, AABB& b, real padding ) const
{
    (void)model;
    b = box;
    b.pad( padding );
    return true;
}

bool Model::Instance::hit( const Model * model, const real3& origin, const real3& direction, const real3& direction_inv, 
                           real solid_angle, real t_min, real t_max, HitInfo& hit_info )
{
    die_assert( sizeof(real) == 4, "real is not float" );
    die_assert( kind == INSTANCE_KIND::MODEL_PTR, "did not get model_ptr instance" );

    //-----------------------------------------------------------
    // Use inverse matrix to transform origin, direction, and direction_inv into target model's space.
    // Then call model's root BVH hit node.
    //-----------------------------------------------------------
    Model  *   t_model         = u.model_ptr;
    Matrix *   M_inv           = &model->matrixes[matrix_inv_i];
    real3      t_origin, t_direction, t_direction_inv;
    M_inv->transform( origin, t_origin );
    M_inv->transform( direction, t_direction );
    M_inv->transform( direction_inv, t_direction_inv );
    mdout << "Model::Instance::hit: M_inv=" << *M_inv << 
                                  " origin="   << origin   << " direction="   << direction   << " direction_inv="   << direction_inv << 
                                  " t_origin=" << t_origin << " t_direction=" << t_direction << " t_direction_inv=" << t_direction_inv << "\n";

    BVH_Node * t_bvh = &t_model->bvh_nodes[t_model->hdr->bvh_root_i];
    if ( !t_bvh->hit( t_model, t_origin, t_direction, t_direction_inv, solid_angle, t_min, t_max, hit_info ) ) return false;

    //-----------------------------------------------------------
    // Use matrix to transform hit_info.p back to global world space.
    // Use transposed inverse matrix to transform hit_info.normal correctly (ask Pete Shirley).
    //-----------------------------------------------------------
    Matrix * M           = &model->matrixes[matrix_i];
    Matrix * M_inv_trans = &model->matrixes[matrix_inv_trans_i];
    real3 p = hit_info.p;
    real3 normal = hit_info.normal;
    M->transform( p, hit_info.p );
    M_inv_trans->transform( normal, hit_info.normal );
    return true;
}

inline bool Model::BVH_Node::bounding_box( const Model * model, Model::AABB& b ) const
{
    (void)model;
    b = box;
    return true;
}

inline bool Model::BVH_Node::hit( const Model * model, const Model::real3& origin, 
                                  const Model::real3& direction, const Model::real3& direction_inv,
                                  Model::real solid_angle, Model::real t_min, Model::real t_max, Model::HitInfo& hit_info ) const
{
    bool r = false;
    uint bvh_i = this - model->bvh_nodes;
    mdout << "Model::BVH_Node::hit: bvh_i=" << bvh_i << "\n";
    if ( box.hit( origin, direction, direction_inv, t_min, t_max ) ) {
        HitInfo left_hit_info;
        HitInfo right_hit_info;
        bool hit_left  = (left_kind == BVH_NODE_KIND::POLYGON)   ? model->polygons[left_i].hit(    model, origin, direction, direction_inv,   
                                                                                                   solid_angle, t_min, t_max, left_hit_info  ) :
                         (left_kind == BVH_NODE_KIND::INSTANCE)  ? model->instances[left_i].hit(   model, origin, direction, direction_inv, 
                                                                                                   solid_angle, t_min, t_max, left_hit_info  ) :
                                                                   model->bvh_nodes[left_i].hit(   model, origin, direction, direction_inv, 
                                                                                                   solid_angle, t_min, t_max, left_hit_info  );
        bool hit_right = (left_i == right_i)                     ? false :   // lone leaf
                         (right_kind == BVH_NODE_KIND::POLYGON)  ? model->polygons[right_i].hit(   model, origin, direction, direction_inv, 
                                                                                                   solid_angle, t_min, t_max, right_hit_info ) :
                         (right_kind == BVH_NODE_KIND::INSTANCE) ? model->instances[right_i].hit(  model, origin, direction, direction_inv, 
                                                                                                   solid_angle, t_min, t_max, right_hit_info ) :
                                                                   model->bvh_nodes[right_i].hit(  model, origin, direction, direction_inv, 
                                                                                                   solid_angle, t_min, t_max, right_hit_info );
        if ( hit_left && (!hit_right || left_hit_info.t < right_hit_info.t) ) {
            hit_info = left_hit_info;
            r = true;
        } else if ( hit_right ) {
            hit_info = right_hit_info;
            r = true;
        }
    }
    mdout << "Model::BVH_Node::hit: bvh_i=" << bvh_i << " return=" << r << "\n";
    return r;
}

#endif
