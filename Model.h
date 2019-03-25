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
//     4) If you want Model to generate a BVH tree for you, add Model::BVH_TREE::BINARY as the third argument
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

    Model( std::string top_file, MIPMAP_FILTER mipmap_filter=MIPMAP_FILTER::NONE, 
                                 BVH_TREE bvh_tree=BVH_TREE::NONE, bool resolve_models=true );
    Model( std::string model_file, bool is_compressed );
    ~Model(); 

    bool write( std::string file_path, bool is_compressed ); 
    bool replace_materials( std::string mtl_file_path );

    static void dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); // utility

    static const uint VERSION = 0xB0BA1f07; // current version 

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    class real4
    {
    public:
        real c[4];
        
        real4( void ) {}
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
        uint64      char_cnt;               // in strings array

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

        uint64      matrix_cnt;             // in matrixes array
        uint64      inst_cnt;               // in instances array

        uint64      bvh_node_cnt;           // in bvh_nodes array
        uint64      bvh_root_i;             // index of root BVH_Node in bvh_nodes array

        uint64      light_cnt;              // in lights array
        real        lighting_scale;         // not sure what this is yet (default: 1.0)
        real3       ambient_intensity;      // ambient light intensity (default: [0.1, 0.1, 0.1]
        real3       background_color;       // default is black, though irrelevant if there's a sky box
        uint        sky_box_tex_i;          // index in textures array of sky box texture
        real        env_map_intensity_scale;// name says it all
        real        opacity_scale;          // not sure what this is for
        real        shadow_caster_count;    // not sure what this is for

        uint64      camera_cnt;             // in cameras array
        uint64      initial_camera_i;       // initial active camera index (default: 0)
        uint64      frame_cnt;              // in frames array 
        uint64      animation_cnt;          // in animations array
        real        animation_speed;        // divide frame time by this to get real time (default: 1.0)
    };

    // TODO: this is a temporary location for these user_defined variables; they will move into Header on next re-gen
        real        tone_white;             // tone mapping white parameter
        real        tone_key;               // tone mapping key parameter

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
        real3           tangent;
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

        // TODO: add value() or color() routine to compute the color using all of these parameters
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

        void   identity(  void );                       // make this the identity matrix
        void   translate( const real3& translation );   // translate this matrix by a real3
        void   scale(     const real3& scaling );       // scale this matrix by a real3
        void   rotate_xz( double radians );             // rotate by radians in xz plane (yaw)
        void   rotate_yz( double radians );             // rotate by radians in yz plane (pitch)
        void   rotate_xy( double radians );             // rotate by radians in xy plane (roll)

        Matrix operator + ( const Matrix& m ) const;    // add two matrices
        Matrix operator - ( const Matrix& m ) const;    // subtract two matrices
        void   multiply(    double s );                 // multiply this matrix by scalar

        real4  row( uint r ) const;                     // returns row r as a vector
        real4  column( uint c ) const;                  // returns column c as a vector
        void   transform( const real4& v, real4& r ) const; // multiply this matrix (lhs) by vector, returning vector r without div by w
        void   transform( const real3& v, real3& r, bool div_by_w=false ) const; // multiply this matrix (lhs) by vector, returning vector r
        void   transform( const Matrix& m, Matrix& r ) const; // multiply this matrix (lhs) by matrix m, returning matrix r
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

    // structs
    char *              mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays
    char *              strings;
    Object *            objects;
    Polygon *           polygons;
    Vertex *            vertexes;
    real3 *             positions;
    real3 *             normals;
    real2 *             texcoords;
    Material *          materials;
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
            
    char * fsc_start;
    char * fsc_end;
    char * fsc;
    char * obj_start;
    char * obj_end;
    char * obj;
    uint   line_num;

    bool load_fsc( std::string fsc_file, std::string dir_name );        // .fscene 
    bool load_obj( std::string obj_file, std::string dir_name );        // .obj
    bool load_mtl( std::string mtl_file, std::string dir_name, bool replacing=false );
    bool load_tex( const char * tex_name, std::string dir_name, Texture *& texture );

    bool open_and_read( std::string file_name, char *& start, char *& end );

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
    bool parse_option_name( std::string& option_name, char *& xxx, char *& xxx_end );
    bool parse_obj_cmd( obj_cmd_t& cmd );
    bool parse_mtl_cmd( mtl_cmd_t& cmd, char *& mtl, char *& mtl_end );
    bool parse_real3( real3& r3, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real2( real2& r2, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real( real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first=false );
    bool parse_int( int& i, char *& xxx, char *& xxx_end );
    bool parse_uint( uint& u, char *& xxx, char *& xxx_end );
    bool parse_bool( bool& b, char *& xxx, char *& xxx_end );
    std::string surrounding_lines( char *& xxx, char *& xxx_end );

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

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); return false; }
#define fsc_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); goto error;   }
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

Model::Model( std::string top_file, Model::MIPMAP_FILTER mipmap_filter, Model::BVH_TREE bvh_tree, bool resolve_models )
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
    hdr->lighting_scale = 1.0;
    hdr->ambient_intensity = real3( 0.1, 0.1, 0.1 );
    hdr->sky_box_tex_i = uint(-1);
    hdr->env_map_intensity_scale = 1.0;
    hdr->opacity_scale = 1.0;
    hdr->initial_camera_i = uint(-1);
    hdr->animation_speed = 1.0;
    hdr->background_color = real3( 0, 0, 0 );

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->obj_cnt     =   128*1024;
    max->poly_cnt    =  1024*1024;
    max->mtl_cnt     =       128;

    max->vtx_cnt     = 4*max->poly_cnt;
    max->pos_cnt     = max->vtx_cnt;
    max->norm_cnt    = max->vtx_cnt;
    max->texcoord_cnt= max->vtx_cnt;
    max->mipmap_filter= mipmap_filter;
    max->tex_cnt     = max->mtl_cnt;
    max->texel_cnt   = max->mtl_cnt * 128*1024;
    max->char_cnt    = max->obj_cnt * 128;
    max->bvh_node_cnt= max->poly_cnt / 2;
    max->matrix_cnt  = 1;
    max->inst_cnt    = 1;
    max->light_cnt   = 1;
    max->camera_cnt  = 1;
    max->frame_cnt   = 1;
    max->animation_cnt = 1;

    //------------------------------------------------------------
    // Allocate arrays
    //------------------------------------------------------------
    strings         = aligned_alloc<char>(     max->char_cnt );
    objects         = aligned_alloc<Object>(   max->obj_cnt );
    polygons        = aligned_alloc<Polygon>(  max->poly_cnt );
    vertexes        = aligned_alloc<Vertex>(   max->vtx_cnt );
    positions       = aligned_alloc<real3>(    max->pos_cnt );
    normals         = aligned_alloc<real3>(    max->norm_cnt );
    texcoords       = aligned_alloc<real2>(    max->texcoord_cnt );
    materials       = aligned_alloc<Material>( max->mtl_cnt );
    textures        = aligned_alloc<Texture>(  max->tex_cnt );
    texels          = aligned_alloc<unsigned char>( max->texel_cnt );
    bvh_nodes       = aligned_alloc<BVH_Node>( max->bvh_node_cnt );
    matrixes        = aligned_alloc<Matrix>(   max->matrix_cnt );
    instances       = aligned_alloc<Instance>( max->inst_cnt );
    lights          = aligned_alloc<Light>(    max->light_cnt );
    cameras         = aligned_alloc<Camera>(   max->camera_cnt );
    frames          = aligned_alloc<Frame>(    max->frame_cnt );
    animations      = aligned_alloc<Animation>(max->animation_cnt );

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
                    uint64( hdr->bvh_node_cnt ) * sizeof( bvh_nodes[0] ) +
                    uint64( hdr->matrix_cnt   ) * sizeof( matrixes[0] ) +
                    uint64( hdr->inst_cnt     ) * sizeof( instances[0] ) +
                    uint64( hdr->light_cnt    ) * sizeof( lights[0] ) +
                    uint64( hdr->camera_cnt   ) * sizeof( cameras[0] ) +
                    uint64( hdr->frame_cnt    ) * sizeof( frames[0] ) +
                    uint64( hdr->animation_cnt) * sizeof( animations[0] );

    is_good = true;
}

Model::Model( std::string model_path, bool is_compressed )
{
    is_good = false;
    mapped_region = nullptr;
    if ( !is_compressed ) {
        read_uncompressed( model_path );
        return;
    }

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
                    error_msg = "could not gzread() file " + model_path + " - gzread() error: " + strerror( errno ); \
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
        error_msg = "hdr->version does not match VERSION=" + std::to_string(VERSION) + ", got " + std::to_string(hdr->version);
        assert( 0 ); \
        return;
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _read( strings,     char,          hdr->char_cnt );
    _read( objects,     Object,        hdr->obj_cnt );
    _read( polygons,    Polygon,       hdr->poly_cnt );
    _read( vertexes,    Vertex,        hdr->vtx_cnt );
    _read( positions,   real3,         hdr->pos_cnt );
    _read( normals,     real3,         hdr->norm_cnt );
    _read( texcoords,   real2,         hdr->texcoord_cnt );
    _read( materials,   Material,      hdr->mtl_cnt );
    _read( textures,    Texture,       hdr->tex_cnt );
    _read( texels,      unsigned char, hdr->texel_cnt );
    _read( bvh_nodes,   BVH_Node,      hdr->bvh_node_cnt );
    _read( matrixes,    Matrix,        hdr->matrix_cnt );
    _read( instances,   Instance,      hdr->inst_cnt );
    _read( lights,      Light,         hdr->light_cnt );
    _read( cameras,     Camera,        hdr->camera_cnt );
    _read( frames,      Frame,         hdr->frame_cnt );
    _read( animations,  Animation,     hdr->animation_cnt );

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

    _write( hdr,         1                  * sizeof(hdr[0]) );
    _write( strings,     hdr->char_cnt      * sizeof(strings[0]) );
    _write( objects,     hdr->obj_cnt       * sizeof(objects[0]) );
    _write( polygons,    hdr->poly_cnt      * sizeof(polygons[0]) );
    _write( vertexes,    hdr->vtx_cnt       * sizeof(vertexes[0]) );
    _write( positions,   hdr->pos_cnt       * sizeof(positions[0]) );
    _write( normals,     hdr->norm_cnt      * sizeof(normals[0]) );
    _write( texcoords,   hdr->texcoord_cnt  * sizeof(texcoords[0]) );
    _write( materials,   hdr->mtl_cnt       * sizeof(materials[0]) );
    _write( textures,    hdr->tex_cnt       * sizeof(textures[0]) );
    _write( texels,      hdr->texel_cnt     * sizeof(texels[0]) );
    _write( bvh_nodes,   hdr->bvh_node_cnt  * sizeof(bvh_nodes[0]) );
    _write( matrixes,    hdr->matrix_cnt    * sizeof(matrixes[0]) );
    _write( instances,   hdr->inst_cnt      * sizeof(instances[0]) );
    _write( lights,      hdr->light_cnt     * sizeof(lights[0]) );
    _write( cameras,     hdr->camera_cnt    * sizeof(cameras[0]) );
    _write( frames,      hdr->frame_cnt     * sizeof(frames[0]) );
    _write( animations,  hdr->animation_cnt * sizeof(animations[0]) );

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
            assert( old_max_cnt != uint(-1) );
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

    _uwrite( hdr,         1                  * sizeof(hdr[0]) );
    _uwrite( strings,     hdr->char_cnt      * sizeof(strings[0]) );
    _uwrite( objects,     hdr->obj_cnt       * sizeof(objects[0]) );
    _uwrite( polygons,    hdr->poly_cnt      * sizeof(polygons[0]) );
    _uwrite( vertexes,    hdr->vtx_cnt       * sizeof(vertexes[0]) );
    _uwrite( positions,   hdr->pos_cnt       * sizeof(positions[0]) );
    _uwrite( normals,     hdr->norm_cnt      * sizeof(normals[0]) );
    _uwrite( texcoords,   hdr->texcoord_cnt  * sizeof(texcoords[0]) );
    _uwrite( materials,   hdr->mtl_cnt       * sizeof(materials[0]) );
    _uwrite( textures,    hdr->tex_cnt       * sizeof(textures[0]) );
    _uwrite( texels,      hdr->texel_cnt     * sizeof(texels[0]) );
    _uwrite( bvh_nodes,   hdr->bvh_node_cnt  * sizeof(bvh_nodes[0]) );
    _uwrite( matrixes,    hdr->matrix_cnt    * sizeof(matrixes[0]) );
    _uwrite( instances,   hdr->inst_cnt      * sizeof(instances[0]) );
    _uwrite( lights,      hdr->light_cnt     * sizeof(lights[0]) );
    _uwrite( cameras,     hdr->camera_cnt    * sizeof(cameras[0]) );
    _uwrite( frames,      hdr->frame_cnt     * sizeof(frames[0]) );
    _uwrite( animations,  hdr->animation_cnt * sizeof(animations[0]) );

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

bool Model::read_uncompressed( std::string model_path )
{
    char * start;
    char * end;
    if ( !open_and_read( model_path, start, end ) ) return false;
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
        rtn_assert( 0, "hdr->version does not match VERSION=" + std::to_string(VERSION) + ", got " + std::to_string(hdr->version) );
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _uread( strings,     char,          hdr->char_cnt );
    _uread( objects,     Object,        hdr->obj_cnt );
    _uread( polygons,    Polygon,       hdr->poly_cnt );
    _uread( vertexes,    Vertex,        hdr->vtx_cnt );
    _uread( positions,   real3,         hdr->pos_cnt );
    _uread( normals,     real3,         hdr->norm_cnt );
    _uread( texcoords,   real2,         hdr->texcoord_cnt );
    _uread( materials,   Material,      hdr->mtl_cnt );
    _uread( textures,    Texture,       hdr->tex_cnt );
    _uread( texels,      unsigned char, hdr->texel_cnt );
    _uread( bvh_nodes,   BVH_Node,      hdr->bvh_node_cnt );
    _uread( matrixes,    Matrix,        hdr->matrix_cnt );
    _uread( instances,   Instance,      hdr->inst_cnt );
    _uread( lights,      Light,         hdr->light_cnt );
    _uread( cameras,     Camera,        hdr->camera_cnt );
    _uread( frames,      Frame,         hdr->frame_cnt );
    _uread( animations,  Animation,     hdr->animation_cnt );

    is_good = true;

    return true;
}

bool Model::load_fsc( std::string fsc_file, std::string dir_name )
{
    tone_key = real(0.2);
    tone_white = real(3.0);
    (void)dir_name;

    //------------------------------------------------------------
    // Map in .fscene file
    //------------------------------------------------------------
    line_num = 1;
    fsc = nullptr;
    uint initial_camera_name_i = uint(-1);
    if ( !open_and_read( fsc_file, fsc_start, fsc_end ) ) goto error;
    fsc = fsc_start;

    //------------------------------------------------------------
    // Parse top dictionary.
    //------------------------------------------------------------
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

                } else if ( strcmp( field, "tone_white" ) == 0 ) {
                    if ( !parse_real( tone_white, fsc, fsc_end, true ) ) goto error;

                } else if ( strcmp( field, "tone_key" ) == 0 ) {
                    if ( !parse_real( tone_key, fsc, fsc_end, true ) ) goto error;

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
        assert( 0 );
        return false;
}

bool Model::load_obj( std::string obj_file, std::string dir_name )
{
    (void)dir_name;

    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    line_num = 1;
    if ( !open_and_read( obj_file, obj_start, obj_end ) ) return false;
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
        assert( 0 );
        return false;
}

bool Model::load_mtl( std::string mtl_file, std::string dir_name, bool replacing )
{
    //------------------------------------------------------------
    // Map in .obj file
    //------------------------------------------------------------
    char * mtl;
    char * mtl_end;
    if ( dir_name != std::string( "" ) ) mtl_file = dir_name + "/" + mtl_file;
    if ( !open_and_read( mtl_file, mtl, mtl_end ) ) return false;

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
                        if ( !load_tex( tex_name, dir_name, texture ) ) return false;
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

bool Model::load_tex( const char * tex_name, std::string dir_name, Model::Texture *& texture )
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
    if ( dir_name != "" ) tex_name_s = dir_name + "/" + tex_name_s;
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

bool Model::open_and_read( std::string file_path, char *& start, char *& end )
{
    const char * fname = file_path.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) std::cout << "open_and_read() error reading " << file_path << ": " << strerror( errno ) << "\n";
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
            if ( ch == '\n' && (xxx == fsc || xxx == obj) ) line_num++;
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
            if ( ch == '\n' && (xxx == fsc || xxx == obj) ) line_num++;
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
            if ( *xxx == '\n' && (xxx == fsc || xxx == obj) ) line_num++;
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

inline bool Model::parse_string_i( uint& s_i, char *& xxx, char *& xxx_end )
{
    std::string s;
    if ( !parse_string( s, xxx, xxx_end ) ) return false;

    uint s_len = s.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    s_i = hdr->char_cnt;
    char * to_s = &strings[hdr->char_cnt];
    hdr->char_cnt += s_len + 1;
    memcpy( to_s, s.c_str(), s_len+1 );
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
    rtn_assert( vld, "unable to parse real in " + std::string( ((xxx == fsc) ? ".fscene" : (xxx == obj) ? ".obj" : ".mtl") ) + " file " + surrounding_lines( xxx, xxx_end ) );
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
    rtn_assert( vld, "unable to parse int in " + std::string( ((xxx == fsc) ? ".fscene" : (xxx == obj) ? ".obj" : ".mtl") ) + " file " + surrounding_lines( xxx, xxx_end ) );
    return true;
}

inline bool Model::parse_uint( uint& u, char *& xxx, char *& xxx_end )
{
    int i;
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
        assert( hdr->inst_cnt == 0 );
        hdr->bvh_root_i = bvh_node( true, 0, hdr->poly_cnt, 1 );
    } else {
        assert( hdr->inst_cnt != 0 && hdr->poly_cnt == 0 );
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

    assert( m >= first && m <= (first+n)  );
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
        assert( node->box.encloses( new_box ) );
        if ( n == 2 ) {
            if ( for_polys ) {
                polygons[first+1].bounding_box( this, new_box );
            } else {
                instances[first+1].bounding_box( this, new_box );
            }
            assert( node->box.encloses( new_box ) );
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
        assert( node->box.encloses( bvh_nodes[left_i].box ) && node->box.encloses( bvh_nodes[right_i].box ) );
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
    Matrix mr;
    mr.identity();
    mr.m[0][0] = c;
    mr.m[0][1] = s;
    mr.m[1][0] = -s;
    mr.m[1][1] = c;
    Matrix r = *this;
    mr.transform( r, *this );
}

void Model::Matrix::rotate_xz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix mr;
    mr.identity();
    mr.m[0][0] = c;
    mr.m[0][2] = s;
    mr.m[2][0] = -s;
    mr.m[2][2] = c;
    Matrix r = *this;
    mr.transform( r, *this );
}

void Model::Matrix::rotate_yz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix mr;
    mr.identity();
    mr.m[1][1] = c;
    mr.m[1][2] = -s;
    mr.m[2][1] = s;
    mr.m[2][2] = c;
    Matrix r = *this;
    mr.transform( r, *this );
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

void Model::Matrix::transform( const real4& v, real4& r ) const
{
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

void Model::Matrix::transform( const Matrix& M2, Matrix& M3 ) const
{
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

inline bool Model::AABB::hit( const Model::real3& origin, const Model::real3& direction, const Model::real3& direction_inv, 
                              Model::real tmin, Model::real tmax ) const 
{
    (void)direction;
    for( uint a = 0; a < 3; a++ ) 
    {
        real dir_inv = direction_inv.c[a];
        real v0 = (min.c[a] - origin.c[a]) * dir_inv;
        real v1 = (max.c[a] - origin.c[a]) * dir_inv;
        tmin = std::max( tmin, std::min( v0, v1 ) );
        tmax = std::min( tmax, std::max( v0, v1 ) );
    }
    bool r = tmax >= std::max( tmin, real(0.0) );
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

                if ( vertexes[0].vn_i != uint(-1) ) {
                    hit_info.normal = n0*alpha + n1*gamma + n2*beta;
                    hit_info.normal.normalize();
                } else {
                    hit_info.normal = normal;
                }
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

                    float opacity = float(mdata[nchan*i + nchan*mx*j+0]) / 255.0;

                    if (float(uniform()) > opacity) {
                        return false;
                    }
                }

                hit_info.model = model;

                return true;
            }
        }
    }
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
    assert( kind == INSTANCE_KIND::MODEL_PTR );

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
    return r;
}

#endif
