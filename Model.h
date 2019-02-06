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
//        Model * model = new Model( "~/models/sanmiguel", "sanmiguel.obj" );   // model dir, .obj file
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//        model->write( "~models/sanmiguel/sanmiguel.model" );   // will write out the self-contained binary model
//                                                               // default is compressed, add "false" as 2nd arg for uncompressed
//
//     2) After that, you can quickly read in the single binary model file using:
//
//        Model * model = new Model( "~models/sanmiguel/sanmiguel.model" );  // add "false" as 2nd arg for uncompressed
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//
//     3) If you want Model to generate mipmap textures, add Model::MIPMAP_FILTER::BOX as a third argument 
//        to the Model() constructor in (2) to get a box filter.  Other filters may be added later.
//        The texture for mip level 0 is the original texture.  The texels for the other mip levels
//        follow immediately with no padding in between.  The original width and height need not
//        be powers-of-2 or equal to each other.  Each mip level has 
//        width=min(1, prev_width>>1) x height=min(1, prev_height>>1) texels.
//        The last mip level will always contain 1x1 texels.  Model does not currently store
//        the number of texels in each level, so you'll need to compute those on the fly as you
//        try to obtain the starting offset for a given level.
//
//     4) If you want Model to generate a BVH tree for you, add Model::BVH_TREE::BINARY as the fourth argument
//        to the Model() constructor in (2) to get a binary BVH tree.  QUAD and OCT trees are not
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
//     5) Optionally generate BVH tree.
//     6) Write to  uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//     7) Read from uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//        The O/S will do the equivalent of an mmap() for each array.
//
// To Do:
//
//     1) Support rest of .obj and .mtl commands and options.
//        - v: we disregard 4+ extra params
//     2) Support the .fbx and .max binary formats.  Currently, you must find an fbx-to-obj conversion program.
//
#ifndef _Model_h
#define _Model_h

#include <cstdint>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <iostream>

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

class Model
{
public:
    typedef uint32_t uint;                  // by default, we use 32-bit indexes for geometry
    typedef uint64_t uint64;                // by default, we use 64-bit indexes for texels
    typedef float    real;

    enum class MIPMAP_FILTER
    {
        NONE,                               // do not generate mipmap levels
        BOX,                                // generate mipmap levels using simple box filter (average)
    };

    enum class BVH_TREE
    {
        NONE,                               // do not generate BVH tree
        BINARY,                             // generate binary BVH tree
    };

    Model( std::string dir_path, std::string obj_file, MIPMAP_FILTER mipmap_filter=MIPMAP_FILTER::NONE, BVH_TREE bvh_tree=BVH_TREE::NONE );
    Model( std::string file_path, bool is_compressed=true );
    ~Model(); 

    bool write( std::string file_path, bool is_compressed=true ); 
    bool replace_materials( std::string mtl_file_path );

    static const uint VERSION = 0xB0BA1f05; // current version is 5

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    class real3
    {
    public:
        real c[3];
        
        real3( void ) {}
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

        real2( void ) {}
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
        uint64      obj_cnt;                // in objects array  
        uint64      poly_cnt;               // in polygons array 
        uint64      vtx_cnt;                // in vertexes array
        uint64      pos_cnt;                // in positions array
        uint64      norm_cnt;               // in normals array
        uint64      texcoord_cnt;           // in texcoords array
        MIPMAP_FILTER mipmap_filter;        // if NONE, there are no mip levels beyond level 0
        uint64      mtl_cnt;                // in materials array
        uint64      tex_cnt;                // in textures array
        uint64      texel_cnt;              // in texels array  (last in file)
        uint64      char_cnt;               // in strings array
        uint64      bvh_node_cnt;           // in bvh_nodes array
        uint64      bvh_root_i;             // index of root bvh_node in bvh_nodes array
    };

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
        uint        poly_i;
        real        t;  
        real        u;
        real        v;
        real        solid_angle;            // solid angle of hit surface (used by mip_level calculation)
        real3       p;
        real3       normal;
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
        bool hit( const real3& origin, const real3& direction, real tmin, real tmax ) const; 
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
        bool hit( const Model * model, const real3& origin, const real3& direction, real t_min, real t_max, HitInfo& hit_info ) const;
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
        real            Ns;                 // specular exponent (focus of specular highlight)  (default: [0,0,0])
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
    };

    class Texture
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        uint            width;              // level 0 width
        uint            height;             // level 0 height
        uint            nchan;              // number of channels (typically 3 or 4)
        uint64          texel_i;            // index into texels array of first texel
        real            bump_multiplier;    // should be for bump maps only but also used with other textures
    };

    class BVH_Node
    {
    public:
        AABB            box;                // bounding box
        bool            left_is_leaf;       // if true, left is a polygon
        bool            right_is_leaf;      // if true, right is a polygon
        uint            left_i;             // index into bvh_nodes or polygons array of left subtree 
        uint            right_i;            // index into bvh_nodes or polygons array of right subtree 

        bool bounding_box( const Model * model, AABB& b ) const;
        bool hit( const Model * model, const real3& origin, const real3& direction, real t_min, real t_max, HitInfo& hit_info ) const;
    };

    // structs
    char *              mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays
    Object *            objects;
    Polygon *           polygons;
    Vertex *            vertexes;
    real3 *             positions;
    real3 *             normals;
    real2 *             texcoords;
    Material *          materials;
    Texture *           textures;
    BVH_Node *          bvh_nodes;
    unsigned char *     texels;
    char *              strings;

    // maps
    std::map<std::string, uint> name_to_mtl_i;
    std::map<std::string, uint> name_to_tex_i;

    // debug flags
    static bool         debug_hit;
    static uint         debug_tex_i;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
// IMPLEMENTATION
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
            
    char *              obj;
    char *              obj_start;
    char *              obj_end;
    uint                line_num;

    bool load_obj( std::string dir_path, std::string obj_file );
    bool load_mtl( std::string dir_path, std::string mtl_file, bool replacing=false );
    bool load_tex( std::string dir_path, char *& tex_name, Texture *& texture );
    bool open_and_read( std::string dir_path, std::string file_name, char *& start, char *& end );

    void skip_whitespace_to_eol( char *& xxx, char *& xxx_end );  // on this line only
    void skip_whitespace( char *& xxx, char *& xxx_end );
    void skip_to_eol( char *& xxx, char *& xxx_end );
    bool eol( char *& xxx, char *& xxx_end );
    bool expect_char( char ch, char *& xxx, char* xxx_end );
    bool expect_cmd( const char * s, char *& xxx, char *& xxx_end );
    bool parse_name( char *& name, char *& xxx, char *& xxx_end );
    bool parse_option_name( std::string& option_name, char *& xxx, char *& xxx_end );
    bool parse_obj_cmd( obj_cmd_t& cmd );
    bool parse_mtl_cmd( mtl_cmd_t& cmd, char *& mtl, char *& mtl_end );
    bool parse_real3( real3& r3, char *& xxx, char *& xxx_end );
    bool parse_real2( real2& r2, char *& xxx, char *& xxx_end );
    bool parse_real( real& r, char *& xxx, char *& xxx_end );
    bool parse_int( int& i, char *& xxx, char *& xxx_end );
    std::string surrounding_lines( char *& xxx, char *& xxx_end );

    // BVH builder
    void bvh_build( BVH_TREE bvh_tree );
    uint bvh_qsplit( uint poly_i, uint n, real pivot, uint axis );
    uint bvh_node( uint poly_i, uint n, uint axis );

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

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); return false; }
#define obj_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); goto error;   }

inline std::istream& operator >> ( std::istream& is, Model::real3& v ) 
{
    is >> v.c[0] >> v.c[1] >> v.c[2];
    return is;
}

inline std::ostream& operator << ( std::ostream& os, const Model::real3& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "," << v.c[2] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::AABB& box ) 
{
    os << box.min << ".." << box.max;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Model::BVH_Node& bvh ) 
{
    os << "box=" << bvh.box << " left_is_leaf=" << bvh.left_is_leaf << " right_is_leaf=" << bvh.right_is_leaf <<
         " left_i=" << bvh.left_i << " right_i=" << bvh.right_i;
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

inline std::ostream& operator << ( std::ostream& os, const Model::Texture& tex ) 
{
    os << "name_i=" << tex.name_i << " width=" << tex.width << " height=" << tex.height << " nchan=" << tex.nchan <<
                    " texel_i=" << tex.texel_i << " bump_multiplier=" << tex.bump_multiplier;
    return os;
}

Model::Model( std::string dir_path, std::string top_file, Model::MIPMAP_FILTER mipmap_filter, Model::BVH_TREE bvh_tree )
{
    is_good = false;
    error_msg = "<unknown error>";
    mapped_region = nullptr;

    hdr = aligned_alloc<Header>( 1 );
    memset( hdr, 0, sizeof( Header ) );
    hdr->version = VERSION;
    hdr->pos_cnt = 1;
    hdr->norm_cnt = 1;
    hdr->texcoord_cnt = 1;
    hdr->mipmap_filter = mipmap_filter;
    hdr->bvh_node_cnt = 0;
    hdr->bvh_root_i = 0;

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->obj_cnt     =   128*1024;
    max->poly_cnt    =  1024*1024;
    max->mtl_cnt     =       1024;

    max->vtx_cnt     = 4*max->poly_cnt;
    max->pos_cnt     = max->vtx_cnt;
    max->norm_cnt    = max->vtx_cnt;
    max->texcoord_cnt= max->vtx_cnt;
    max->mipmap_filter= mipmap_filter;
    max->tex_cnt     = max->mtl_cnt;
    max->texel_cnt   = max->mtl_cnt * 128*1024;
    max->char_cnt    = max->obj_cnt * 128;
    max->bvh_node_cnt= max->poly_cnt / 2;

    //------------------------------------------------------------
    // Allocate arrays
    //------------------------------------------------------------
    objects         = aligned_alloc<Object>(   max->obj_cnt );
    polygons        = aligned_alloc<Polygon>(  max->poly_cnt );
    vertexes        = aligned_alloc<Vertex>(   max->vtx_cnt );
    positions       = aligned_alloc<real3>(    max->pos_cnt );
    normals         = aligned_alloc<real3>(    max->norm_cnt );
    texcoords       = aligned_alloc<real2>(    max->texcoord_cnt );
    materials       = aligned_alloc<Material>( max->mtl_cnt );
    textures        = aligned_alloc<Texture>(  max->tex_cnt );
    texels          = aligned_alloc<unsigned char>( max->texel_cnt );
    strings         = aligned_alloc<char>(     max->char_cnt );
    bvh_nodes       = aligned_alloc<BVH_Node>( max->bvh_node_cnt );

    //------------------------------------------------------------
    // Load .json or .obj depending on file extension
    //------------------------------------------------------------
    if ( !load_obj( dir_path, top_file ) ) return;
    is_good = true;

    //------------------------------------------------------------
    // Optionally build BVH
    //------------------------------------------------------------
    if ( bvh_tree != BVH_TREE::NONE ) bvh_build( bvh_tree );
}

Model::Model( std::string file_path, bool is_compressed )
{
    is_good = false;
    mapped_region = nullptr;
    if ( !is_compressed ) {
        read_uncompressed( file_path );
        return;
    }

    gzFile fd = gzopen( file_path.c_str(), "r" );
    if ( fd == Z_NULL ) {
        "Could not gzopen() file " + file_path + " for reading - gzopen() error: " + strerror( errno );
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
                assert( 0 ); \
                return; \
            } \
            char * _addr = reinterpret_cast<char *>( array ); \
            for( uint64 _byte_cnt = (cnt)*sizeof(type); _byte_cnt != 0;  ) \
            { \
                uint64 _this_byte_cnt = 1024*1024*1024; \
                if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
                if ( gzread( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                    gzclose( fd ); \
                    error_msg = "could not gzread() file " + file_path + " - gzread() error: " + strerror( errno ); \
                    assert( 0 ); \
                    return; \
                } \
                _byte_cnt -= _this_byte_cnt; \
                _addr     += _this_byte_cnt; \
            } \
        } \

    _read( hdr,         Header,   1 );
    if ( hdr->version != VERSION ) {
        gzclose( fd );
        error_msg = "hdr->version does not match VERSION";
        assert( 0 ); \
        return;
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _read( objects,     Object,   hdr->obj_cnt );
    _read( polygons,    Polygon,  hdr->poly_cnt );
    _read( vertexes,    Vertex,   hdr->vtx_cnt );
    _read( positions,   real3,    hdr->pos_cnt );
    _read( normals,     real3,    hdr->norm_cnt );
    _read( texcoords,   real2,    hdr->texcoord_cnt );
    _read( materials,   Material, hdr->mtl_cnt );
    _read( textures,    Texture,  hdr->tex_cnt );
    _read( texels,      unsigned char, hdr->texel_cnt );
    _read( strings,     char,     hdr->char_cnt );
    _read( bvh_nodes,   BVH_Node, hdr->bvh_node_cnt );

    gzclose( fd );

    is_good = true;
}

Model::~Model()
{
    if ( mapped_region != nullptr ) {
        delete mapped_region;
        mapped_region = nullptr;
    } else {
        delete objects;
        delete polygons;
        delete vertexes;
        delete positions;
        delete normals;
        delete texcoords;
        delete materials;
        delete textures;
        delete texels;
        delete strings;
        delete bvh_nodes;
    }
}

bool Model::write( std::string file_path, bool is_compressed ) 
{
    if ( !is_compressed ) return write_uncompressed( file_path );

    gzFile fd = gzopen( file_path.c_str(), "w" );
    rtn_assert( fd != Z_NULL, "could not gzopen() file " + file_path + " for writing - gzopen() error: " + strerror( errno ) );

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
                rtn_assert( 0, "could not gzwrite() file " + file_path + " - gzwrite() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _write( hdr,         1                 * sizeof(hdr[0]) );
    _write( objects,     hdr->obj_cnt      * sizeof(objects[0]) );
    _write( polygons,    hdr->poly_cnt     * sizeof(polygons[0]) );
    _write( vertexes,    hdr->vtx_cnt      * sizeof(vertexes[0]) );
    _write( positions,   hdr->pos_cnt      * sizeof(positions[0]) );
    _write( normals,     hdr->norm_cnt     * sizeof(normals[0]) );
    _write( texcoords,   hdr->texcoord_cnt * sizeof(texcoords[0]) );
    _write( materials,   hdr->mtl_cnt      * sizeof(materials[0]) );
    _write( textures,    hdr->tex_cnt      * sizeof(textures[0]) );
    _write( texels,      hdr->texel_cnt    * sizeof(texels[0]) );
    _write( strings,     hdr->char_cnt     * sizeof(strings[0]) );
    _write( bvh_nodes,   hdr->bvh_node_cnt * sizeof(bvh_nodes[0]) );

    gzclose( fd );
    return true;
}

bool Model::replace_materials( std::string mtl_file_path )
{
    return load_mtl( "", mtl_file_path, true );
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
            assert( old_max_cnt != uint(-1) );
            max_cnt = uint(-1);
        }
        posix_memalign( &mem, getpagesize(), max_cnt*sizeof(T) );
        memcpy( mem, array, hdr_cnt*sizeof(T) );
        delete array;
        array = reinterpret_cast<T *>( mem );
    }
}

bool Model::write_uncompressed( std::string file_path ) 
{
    int fd = open( file_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( fd < 0 ) std::cout << "open() for write error: " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open() file " + file_path + " for writing - open() error: " + strerror( errno ) );

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
                rtn_assert( 0, "could not write() file " + file_path + " - write() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _uwrite( hdr,         1                 * sizeof(hdr[0]) );
    _uwrite( objects,     hdr->obj_cnt      * sizeof(objects[0]) );
    _uwrite( polygons,    hdr->poly_cnt     * sizeof(polygons[0]) );
    _uwrite( vertexes,    hdr->vtx_cnt      * sizeof(vertexes[0]) );
    _uwrite( positions,   hdr->pos_cnt      * sizeof(positions[0]) );
    _uwrite( normals,     hdr->norm_cnt     * sizeof(normals[0]) );
    _uwrite( texcoords,   hdr->texcoord_cnt * sizeof(texcoords[0]) );
    _uwrite( materials,   hdr->mtl_cnt      * sizeof(materials[0]) );
    _uwrite( textures,    hdr->tex_cnt      * sizeof(textures[0]) );
    _uwrite( texels,      hdr->texel_cnt    * sizeof(texels[0]) );
    _uwrite( strings,     hdr->char_cnt     * sizeof(strings[0]) );
    _uwrite( bvh_nodes,   hdr->bvh_node_cnt * sizeof(bvh_nodes[0]) );

    fsync( fd ); // flush
    close( fd );

    if ( 0 ) {
        for( uint32_t i = 0; i < hdr->bvh_node_cnt; i++ )
        {
            std::cout << "bvh[" << i << "]: " << bvh_nodes[i] << "\n";
        }
    }
    return true;
}

bool Model::read_uncompressed( std::string file_path )
{
    char * start;
    char * end;
    if ( !open_and_read( "", file_path, start, end ) ) return false;
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
        rtn_assert( 0, "hdr->version does not match VERSION" );
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _uread( objects,     Object,   hdr->obj_cnt );
    _uread( polygons,    Polygon,  hdr->poly_cnt );
    _uread( vertexes,    Vertex,   hdr->vtx_cnt );
    _uread( positions,   real3,    hdr->pos_cnt );
    _uread( normals,     real3,    hdr->norm_cnt );
    _uread( texcoords,   real2,    hdr->texcoord_cnt );
    _uread( materials,   Material, hdr->mtl_cnt );
    _uread( textures,    Texture,  hdr->tex_cnt );
    _uread( texels,      unsigned char, hdr->texel_cnt );
    _uread( strings,     char,     hdr->char_cnt );
    _uread( bvh_nodes,   BVH_Node, hdr->bvh_node_cnt );

    is_good = true;

    return true;
}

bool Model::load_obj( std::string dir_path, std::string obj_file )
{
    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    if ( !open_and_read( dir_path, obj_file, obj_start, obj_end ) ) return false;
    obj = obj_start;
    line_num = 1;

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
        if ( obj == obj_end ) {
            // done, no errors
            hdr->byte_cnt = uint64( 1                 ) * sizeof( hdr ) +
                            uint64( hdr->obj_cnt      ) * sizeof( objects[0] ) +
                            uint64( hdr->poly_cnt     ) * sizeof( polygons[0] ) +
                            uint64( hdr->vtx_cnt      ) * sizeof( vertexes[0] ) +
                            uint64( hdr->pos_cnt      ) * sizeof( positions[0] ) +
                            uint64( hdr->norm_cnt     ) * sizeof( normals[0] ) +
                            uint64( hdr->texcoord_cnt ) * sizeof( texcoords[0] ) +
                            uint64( hdr->mtl_cnt      ) * sizeof( materials[0] ) +
                            uint64( hdr->tex_cnt      ) * sizeof( textures[0] ) +
                            uint64( hdr->texel_cnt    ) * sizeof( texels[0] ) +
                            uint64( hdr->char_cnt     ) * sizeof( strings[0] ) + 
                            uint64( hdr->bvh_node_cnt ) * sizeof( bvh_nodes[0] );
            return true;
        }

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
                polygon = &polygons[ hdr->poly_cnt++ ];
                polygon->mtl_i = mtl_i;
                polygon->vtx_cnt = 0;
                polygon->vtx_i = hdr->vtx_cnt;
                while( !eol( obj, obj_end ) ) 
                {
                    polygon->vtx_cnt++;
                    perhaps_realloc<Vertex>( vertexes, hdr->vtx_cnt, max->vtx_cnt, 1 );
                    vertex = &vertexes[ hdr->vtx_cnt++ ];

                    int v_i;
                    if ( !parse_int( v_i, obj, obj_end ) )   goto error;
                    dprint( "v_i=" + std::to_string( v_i ) );
                    vertex->v_i = (v_i >= 0)  ? v_i : (hdr->pos_cnt + v_i);

                    if ( obj != obj_end && *obj == '/' ) {
                        if ( !expect_char( '/', obj, obj_end ) ) goto error;
                        int vt_i;
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
                            int vn_i;
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
                }
                break;
            }

            case CMD_MTLLIB:
                if ( !parse_name( mtllib, obj, obj_end ) ) goto error;
                if ( !load_mtl( dir_path, mtllib ) ) goto error;
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
        assert( 0 );
        return false;
}

bool Model::load_mtl( std::string dir_path, std::string mtl_file, bool replacing )
{
    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    char * mtl;
    char * mtl_end;
    if ( !open_and_read( dir_path, mtl_file, mtl, mtl_end ) ) return false;

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
                    while( parse_option_name( option_name, mtl, mtl_end ) )
                    {
                        if ( option_name == std::string("-bm") ) {
                            if ( !parse_real( bump_multiplier, mtl, mtl_end ) ) return false;
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
                        if ( !load_tex( dir_path, tex_name, texture ) ) return false;
                    }
                    texture->bump_multiplier = bump_multiplier;

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

bool Model::load_tex( std::string dir_path, char *& tex_name, Model::Texture *& texture )
{
    // allocate Texture structure
    //
    perhaps_realloc<Texture>( textures, hdr->tex_cnt, max->tex_cnt, 1 );
    uint tex_i = hdr->tex_cnt++;
    texture = &textures[tex_i];
    memset( texture, 0, sizeof( Texture ) );
    texture->name_i = tex_name - strings;
    name_to_tex_i[ tex_name ] = tex_i;

    // change any '\' to '/'
    //
    std::string tex_name_s = std::string( tex_name );
    if ( dir_path != "" ) tex_name_s = dir_path + "/" + tex_name_s;
    char * file_name = strdup( tex_name_s.c_str() );
    size_t len = strlen( file_name );
    for( size_t i = 0; i < len; i++ ) 
    {
        if ( file_name[i] == '\\' ) file_name[i] = '/';
    }

    // read file using stb_image.h
    //
    int    w;
    int    h;
    int    nchan;
    unsigned char * data = stbi_load( file_name, &w, &h, &nchan, 0 );
    rtn_assert( data != nullptr, "unable to read in texture file " + std::string( file_name ) );
    texture->width  = w;
    texture->height = h;
    texture->nchan  = nchan;
    texture->texel_i = hdr->texel_cnt;
    uint width      = texture->width;
    uint height     = texture->height;
    uint byte_width = width * texture->nchan;
    uint byte_cnt   = byte_width * height;
    perhaps_realloc<unsigned char>( texels, hdr->texel_cnt, max->texel_cnt, byte_cnt );
    hdr->texel_cnt += byte_cnt;
    memcpy( texels + texture->texel_i, data, byte_cnt );
    delete data;

    // if requested, generate mipmap levels down to 1x1
    //
    uint prior_byte_cnt = 0;
    while( hdr->mipmap_filter != MIPMAP_FILTER::NONE && !(width == 1 && height == 1) )
    {
        uint to_width  = width  >> 1;
        uint to_height = height >> 1;
        if ( to_width  == 0 ) to_width  = 1;
        if ( to_height == 0 ) to_height = 1;

        uint   to_byte_width = to_width * texture->nchan;
        uint   to_byte_cnt   = to_byte_width * to_height;
        perhaps_realloc<unsigned char>( texels, hdr->texel_cnt, max->texel_cnt, to_byte_cnt );
        hdr->texel_cnt += to_byte_cnt;

        unsigned char * rgb  = texels + texture->texel_i + prior_byte_cnt;  // must do this after realloc
        unsigned char * trgb = rgb + byte_cnt;

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
                        unsigned char * frgb = rgb + fii*byte_width + fjj*texture->nchan;
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

        prior_byte_cnt += byte_cnt;
        height          = to_height;
        width           = to_width;
        byte_width      = to_byte_width;
        byte_cnt        = to_byte_cnt;
    }

    return true;
}

bool Model::open_and_read( std::string dir_path, std::string file_name, char *& start, char *& end )
{
    std::string file_path = (dir_path != "") ? (dir_path + "/" + file_name) : file_name;
    const char * fname = file_path.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) std::cout << "open_and_read() error reading " << dir_path << "/" << file_name << ": " << strerror( errno ) << "\n";
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
            rtn_assert( 0, "could not read() file " + std::string(fname) + " - read error: " + std::string( strerror( errno ) ) );
        }
        size -= _this_size;
        addr += _this_size;
    }
    close( fd );
    return true;
}

inline void Model::skip_whitespace( char *& xxx, char *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) return;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx == obj ) line_num++;
            in_comment = false;
        }
        xxx++;
    }
}

inline void Model::skip_whitespace_to_eol( char *& xxx, char *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) return;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx == obj ) line_num++;
            break;
        }
        xxx++;
    }
}

inline void Model::skip_to_eol( char *& xxx, char *& xxx_end )
{
    if ( !eol( xxx, xxx_end ) ) {
        while( xxx != xxx_end )
        {
            char ch = *xxx;
            if ( ch == '\n' || ch == '\r' ) break;
            xxx++;
        }
    }
}

inline bool Model::eol( char *& xxx, char *& xxx_end )
{
    skip_whitespace_to_eol( xxx, xxx_end );

    if ( xxx == xxx_end || *xxx == '\n' || *xxx == '\r' ) {
        if ( xxx != xxx_end ) {
            if ( *xxx == '\n' && xxx == obj ) line_num++;
            xxx++;
        }
        dprint( "at eol" );
        return true;
    } else {
        dprint( "not at eol, char='" + std::string( 1, *xxx ) + "'" );
        return false;
    }
}

inline bool Model::expect_char( char ch, char *& xxx, char* xxx_end )
{
    rtn_assert( xxx != xxx_end, "premature end of file" );
    rtn_assert( *xxx == ch, "unexpected character: " + std::string( 1, *xxx ) );
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

inline bool Model::parse_real3( Model::real3& r3, char *& xxx, char *& xxx_end )
{
    return parse_real( r3.c[0], xxx, xxx_end ) && 
           parse_real( r3.c[1], xxx, xxx_end ) && 
           parse_real( r3.c[2], xxx, xxx_end );
}

inline bool Model::parse_real2( Model::real2& r2, char *& xxx, char *& xxx_end )
{
    return parse_real( r2.c[0], xxx, xxx_end ) && 
           parse_real( r2.c[1], xxx, xxx_end );
}

inline bool Model::parse_real( Model::real& r, char *& xxx, char *& xxx_end )
{
    bool vld = false;
    bool is_neg = false;
    bool in_frac = false;
    bool has_exp = false;
    uint u = 0;     // integer part
    uint f = 0;     // frac part before divide
    int e10 = 0;
    double f_factor = 1.0;
    dprint( "parse_real *xxx=" + std::string( 1, *xxx ) );
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
            is_neg = true;
            xxx++;
            continue;
        }

        if ( ch == '.' && !in_frac ) {
            in_frac = true;
            xxx++;
            continue;
        }

        if ( ch == 'e' || ch == 'E' ) {
            rtn_assert( !has_exp, "real has more than one 'e' exponent" );
            has_exp = true;
            xxx++;
            if ( !parse_int( e10, xxx, xxx_end ) ) return false;
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;

        uint digit = ch - '0';
        if ( in_frac ) {
            f = 10*f + digit;
            f_factor *= 10.0;
        } else {
            u = 10*u + digit;
        }

        vld = true;
        xxx++;
    }

    r = double(u) + f/f_factor;
    if ( is_neg ) r = -r;
    if ( e10 != 0 ) r *= pow( 10.0, e10 );
    dprint( "real=" + std::to_string( r ) );
    rtn_assert( vld, "unable to parse real in " + std::string( ((xxx == obj) ? ".obj" : ".mtl") ) + " file " + surrounding_lines( xxx, xxx_end ) );
    return vld;
}

inline bool Model::parse_int( int& i, char *& xxx, char *& xxx_end )
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
    rtn_assert( vld, "unable to parse int" );
    return true;
}

std::string Model::surrounding_lines( char *& xxx, char *& xxx_end )
{
    uint eol_cnt = 0;
    std::string s = "";
    while( eol_cnt != 2 && xxx != xxx_end )
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
    hdr->bvh_root_i = bvh_node( 0, hdr->poly_cnt, 1 );
}

inline Model::uint Model::bvh_qsplit( Model::uint poly_i, Model::uint n, Model::real pivot, Model::uint axis )
{
   uint m = poly_i;

   for( uint i = poly_i; i < (poly_i+n); i++ )
   {
       AABB box;
       polygons[i].bounding_box( this, box );
       real centroid = (box.min.c[axis] + box.max.c[axis]) * 0.5;
       if ( centroid < pivot ) {
           Polygon temp = polygons[i];
           polygons[i]  = polygons[m];
           polygons[m]  = temp;
           m++;
       }
    }

    assert( m >= poly_i && m <= (poly_i+n)  );
    if ( m == poly_i || m == (poly_i+n) ) m = poly_i + n/2;
    return m;
}

Model::uint Model::bvh_node( Model::uint poly_i, Model::uint n, Model::uint axis ) 
{
    perhaps_realloc( bvh_nodes, hdr->bvh_node_cnt, max->bvh_node_cnt, 1 );
    uint bvh_i = hdr->bvh_node_cnt++;
    BVH_Node * node = &bvh_nodes[bvh_i];

    polygons[poly_i].bounding_box( this, node->box );
    for( uint i = 1; i < n; i++ )
    {
        AABB new_box;
        polygons[poly_i+i].bounding_box( this, new_box );
        node->box.expand( new_box );
    }

    if ( n == 1 || n == 2 ) {
        node->left_is_leaf = true;
        node->right_is_leaf = true;
        node->left_i  = poly_i;
        AABB new_box;
        polygons[poly_i+0].bounding_box( this, new_box );
        assert( node->box.encloses( new_box ) );
        if ( n == 2 ) {
            polygons[poly_i+1].bounding_box( this, new_box );
            assert( node->box.encloses( new_box ) );
            node->right_i = poly_i + 1;
        } else {
            node->right_i = poly_i;
        }

    } else {
        node->left_is_leaf = false;
        node->right_is_leaf = false;
        real pivot = (node->box.min.c[axis] + node->box.max.c[axis]) * 0.5;
        uint m = bvh_qsplit( poly_i, n, pivot, axis );
        uint nm = m - poly_i;
        uint left_i  = bvh_node( poly_i, nm,   (axis + 1) % 3 );
        uint right_i = bvh_node(      m, n-nm, (axis + 1) % 3 );
        node = &bvh_nodes[bvh_i];  // could change after previous calls
        node->left_i  = left_i;
        node->right_i = right_i;
        assert( node->box.encloses( bvh_nodes[left_i].box ) && node->box.encloses( bvh_nodes[right_i].box ) );
    }

    return bvh_i;
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

bool Model::debug_hit = false;

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

inline bool Model::AABB::hit( const Model::real3& origin, const Model::real3& direction, Model::real tmin, Model::real tmax ) const 
{
    for( uint a = 0; a < 3; a++ ) 
    {
        real dir = direction.c[a];
        real v0 = (min.c[a] - origin.c[a]) / dir;
        real v1 = (max.c[a] - origin.c[a]) / dir;
        tmin = std::max( tmin, std::min( v0, v1 ) );
        tmax = std::min( tmax, std::max( v0, v1 ) );
    }
    bool r = tmax >= std::max( tmin, real(0.0) );
    return r;
}

inline bool Model::Polygon::hit( const Model * model, const real3& origin, const real3& direction, real t_min, real t_max, 
                                 HitInfo& hit_info ) const
{
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
                hit_info.poly_i = this - model->polygons;
                hit_info.t = t;
                hit_info.p = p;

                if ( vertexes[0].vn_i != uint(-1) ) {
                    hit_info.normal = n0*alpha + n1*gamma + n2*beta;
                    hit_info.normal.normalize();
                } else {
                    hit_info.normal = normal;
                }

                real distance_squared = t*t / direction.length_sqr();
                hit_info.solid_angle  = area / distance_squared;                                       // off by 2X, but...
                hit_info.solid_angle *= std::abs(0.5*(u0*v1 + u1*v2 + u2*v0 - u0*v2 - u1*v0 - u2*v1)); // correction for UV at corners as lengths
                hit_info.u = alpha*u0 + gamma*u1 + beta*u2 ;
                hit_info.v = alpha*v0 + gamma*v1 + beta*v2 ;
                return true;
             }
        }
    }
    return false;
}

inline bool Model::BVH_Node::bounding_box( const Model * model, Model::AABB& b ) const
{
    (void)model;
    b = box;
    return true;
}

inline bool Model::BVH_Node::hit( const Model * model, const Model::real3& origin, const Model::real3& direction, 
                                  Model::real t_min, Model::real t_max, Model::HitInfo& hit_info ) const
{
    bool r = false;
    if ( box.hit( origin, direction, t_min, t_max ) ) {
        HitInfo left_hit_info;
        HitInfo right_hit_info;
        bool hit_left  = left_is_leaf  ? model->polygons[left_i].hit(   model, origin, direction, t_min, t_max, left_hit_info ) 
                                       : model->bvh_nodes[left_i].hit(  model, origin, direction, t_min, t_max, left_hit_info );
        bool hit_right = (left_i == right_i) ? false :  // lone leaf
                         right_is_leaf ? model->polygons[right_i].hit(  model, origin, direction, t_min, t_max, right_hit_info )
                                       : model->bvh_nodes[right_i].hit( model, origin, direction, t_min, t_max, right_hit_info );
        if ( hit_left && (!hit_right || left_hit_info.t < right_hit_info.t) ) {
            hit_info = left_hit_info;
            r = true;
        } else if ( hit_right ) {
            hit_info = right_hit_info;
            r = true;
        }
    }
    return r;
}

#endif
