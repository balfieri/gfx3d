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
//     1) Please first convert .jpg/.png/etc. textures to .bmp format: 
//
//              cd textures
//              mogrify -format bmp *.{jpg,png,TGA} 
//
//        This way, we don't create a dependency on other code.
//        Model will convert .jpg/.png/etc. suffixes to .bmp automatically, so no need to change .mtl files.
//
//     2) #include "Model.h"
//
//        Model * model = new Model( "~/models/sanmiguel", "sanmiguel.obj" );   // model dir, .obj file
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//        model->write( "~models/sanmiguel/sanmiguel.model" );   // will write out the self-contained binary model
//                                                               // default is compressed, add "false" as 2nd arg for uncompressed
//
//     3) After that, you can quickly read in the single binary model file using:
//
//        Model * model = new Model( "~models/sanmiguel/sanmiguel.model" );  // add "false" as 2nd arg for uncompressed
//        if ( !model->is_good ) {
//            std::cout << "Model load failed with error: " << model->error_msg << "\n";
//            exit( 1 );
//        }
//
//     4) If you want Model to generate mipmap textures, add Model::MIPMAP_FILTER::BOX as a third argument 
//        to the Model() constructor in (2) to get a box filter.  Other filters may be added later.
//        The texture for mip level 0 is the original texture.  The texels for the other mip levels
//        follow immediately with no padding in between.  The original width and height need not
//        be powers-of-2 or equal to each other.  Each mip level has 
//        width=min(1, prev_width>>1) x height=min(1, prev_height>>1) texels.
//        The last mip level will always contain 1x1 texels.  Model does not currently store
//        the number of texels in each level, so you'll need to compute those on the fly as you
//        try to obtain the starting offset for a given level.
//
//     5) If you want Model to generate a BVH tree for you, add Model::BVH_TREE::BINARY as the fourth argument
//        to the Model() constructor in (2) to get a binary BVH tree.  QUAD and OCT trees are not
//        currently supported
//
// How it works:
//
//     1) Preallocate large virtual memory 1D arrays for materials, texels, positions, normals, vertexes, polygons
//     2) Read entire .obj file into memory (the o/s should effectively make this work like an mmap).
//     3) Parse .obj file using custom parser that goes character-by-character and does its own number conversions. 
//     4) Add elements to 1D arrays.  Load any .mtl file encounted in .obj file.
//
// To Do:
//
//     1) Support rest of .obj and .mtl commands and options.
//     2) Support the .fbx binary format.
//     3) This code is set up to write out everything to one binary file.
//        There are no pointers in the data, only indexes.
//        The best way to do it is to put everything in one mmap() region that is backed by a file.
//        Then later mmap() in the same file to get back the same data without even touching it.
//
#ifndef _Model_h
#define _Model_h

#include <cstdint>
#include <string>
#include <cmath>
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

class Model
{
public:
    typedef uint32_t uint;
    typedef float    real;

    class real3
    {
    public:
        real3( void ) {}
        real3( real c0, real c1, real c2 ) { c[0] = c0; c[1] = c1; c[2] = c2; }

        real            c[3];
    };

    class real2
    {
    public:
        real            c[2];
    };

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

    class Header                            // header (of future binary file)
    {
    public:
        uint        version;                // version
        uint64_t    byte_cnt;               // total in-memory bytes including this header
        uint        obj_cnt;                // in objects array  
        uint        poly_cnt;               // in polygons array 
        uint        vtx_cnt;                // in vertexes array
        uint        pos_cnt;                // in positions array
        uint        norm_cnt;               // in normals array
        uint        texcoord_cnt;           // in texcoords array
        MIPMAP_FILTER mipmap_filter;        // if NONE, there are no mip levels beyond level 0
        uint        mtl_cnt;                // in materials array
        uint        tex_cnt;                // in textures array
        uint        texel_cnt;              // in texels array  (last in file)
        uint        char_cnt;               // in strings array
        uint        bvh_node_cnt;           // in bvh_nodes array
        uint        bvh_root_i;             // index of root bvh_node in bvh_nodes array
    };

    class Object
    {
    public:
        uint        name_i;                 // index of object name in strings array
        uint        poly_cnt;               // number of polygons in this object
        uint        poly_i;                 // index of first polygon in polygons array
    };

    class Polygon
    {
    public:
        uint        mtl_i;                  // index into materials array
        uint        vtx_cnt;                // number of vertices
        uint        vtx_i;                  // index into vertexes array of first vertex
        real3       normal;                 // surface normal
        real        area;                   // surface area
    };

    class Vertex
    {
    public:
        uint        v_i;                    // index into positions array
        uint        vn_i;                   // index into normals array
        uint        vt_i;                   // index into texcoords array
    };

    class Material
    {
    public:
        uint            name_i;             // index into strings array
        real3           Ka;
        real3           Kd;
        real3           Ke;
        real3           Ks;
        real3           Tf;
        real            Tr;
        real            Ns;
        real            Ni;
        real            d;
        real            illum;
        uint            map_Ka_i;           // index into textures array
        uint            map_Kd_i;
        uint            map_Ke_i;
        uint            map_Ks_i;
        uint            map_Bump_i;
    };

    class Texture
    {
    public:
        uint            name_i;             // index in strings array (null-terminated strings)
        uint            width;              // level 0 width
        uint            height;             // level 0 height
        uint            texel_i;            // index into texels array of first texel
    };

    class AABB                              // axis aligned bounding box
    {
    public:
        real3           min;                // bounding box min
        real3           max;                // bounding box max

        void expand( const AABB& other )
        {
            for( uint i = 0; i < 3; i++ )
            {
                if ( other.min.c[i] < min.c[i] ) min.c[i] = other.min.c[i];
                if ( other.max.c[i] > max.c[i] ) max.c[i] = other.max.c[i];
            }
        }
    };

    class BVH_Node
    {
    public:
        AABB            box;                // bounding box
        bool            left_is_leaf;       // if true, left is a polygon
        bool            right_is_leaf;      // if true, right is a polygon
        uint            left_i;             // index into bvh_nodes or polygons array of left subtree 
        uint            right_i;            // index into bvh_nodes or polygons array of right subtree 
    };

    // public fields
    //
    static const uint VERSION = 0xB0BA1f03; // current version is 3

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

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
    char *              texels;
    char *              strings;

    // maps
    std::map<std::string, Material *> name_to_mtl;
    std::map<std::string, Texture  *> name_to_tex;



    //------------------------------------------------------------------------------
    // IMPLEMENTATION
    //------------------------------------------------------------------------------

    #define dprint( msg )
    //#define dprint( msg ) std::cout << (msg) << "\n"

    // these are done as macros to avoid evaluating msg (it makes a big difference)
    #define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = msg; assert( false ); return false; }
    #define obj_assert( bool, msg ) if ( !(bool) ) { error_msg = msg; assert( false ); goto error;   }

    // returns array of T on a page boundary
    template<typename T>
    T * aligned_alloc( uint cnt )
    {
        void * mem = nullptr;
        posix_memalign( &mem, getpagesize(), cnt*sizeof(T) );
        return reinterpret_cast<T *>( mem );
    }

    // reallocate array if we are about to exceed its current size
    template<typename T>
    inline void perhaps_realloc( T *& array, uint& hdr_cnt, uint& max_cnt, uint add_cnt )
    {
        while( (hdr_cnt + add_cnt) > max_cnt ) {
            void * mem = nullptr;
            uint old_max_cnt = max_cnt;
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

    Model( std::string dir_path, std::string obj_file, MIPMAP_FILTER mipmap_filter=MIPMAP_FILTER::NONE, BVH_TREE bvh_tree=BVH_TREE::NONE )
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
        texels          = aligned_alloc<char>(     max->texel_cnt );
        strings         = aligned_alloc<char>(     max->char_cnt );
        bvh_nodes       = aligned_alloc<BVH_Node>( max->bvh_node_cnt );

        //------------------------------------------------------------
        // Map in .obj file
        //------------------------------------------------------------
        if ( !open_and_read( dir_path, obj_file, obj_start, obj_end ) ) return;
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
                hdr->byte_cnt = uint64_t( 1                 ) * sizeof( hdr ) +
                                uint64_t( hdr->obj_cnt      ) * sizeof( objects[0] ) +
                                uint64_t( hdr->poly_cnt     ) * sizeof( polygons[0] ) +
                                uint64_t( hdr->vtx_cnt      ) * sizeof( vertexes[0] ) +
                                uint64_t( hdr->pos_cnt      ) * sizeof( positions[0] ) +
                                uint64_t( hdr->norm_cnt     ) * sizeof( normals[0] ) +
                                uint64_t( hdr->texcoord_cnt ) * sizeof( texcoords[0] ) +
                                uint64_t( hdr->mtl_cnt      ) * sizeof( materials[0] ) +
                                uint64_t( hdr->tex_cnt      ) * sizeof( textures[0] ) +
                                uint64_t( hdr->texel_cnt    ) * sizeof( texels[0] ) +
                                uint64_t( hdr->char_cnt     ) * sizeof( strings[0] ) + 
                                uint64_t( hdr->bvh_node_cnt ) * sizeof( bvh_nodes[0] );
                is_good = true;
                if ( bvh_tree != BVH_TREE::NONE ) bvh_build( bvh_tree );
                return;
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
                    
                case CMD_V:
                    perhaps_realloc<real3>( positions, hdr->pos_cnt, max->pos_cnt, 1 );
                    if ( !parse_real3( positions[ hdr->pos_cnt++ ], obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_VN:
                    perhaps_realloc<real3>( normals, hdr->norm_cnt, max->norm_cnt, 1 );
                    if ( !parse_real3( normals[ hdr->norm_cnt++ ], obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_VT:
                    perhaps_realloc<real2>( texcoords, hdr->texcoord_cnt, max->texcoord_cnt, 1 );
                    if ( !parse_real2( texcoords[ hdr->texcoord_cnt++ ], obj, obj_end ) ) goto error;
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

                        if ( !expect_char( '/', obj, obj_end ) ) goto error;

                        int vt_i;
                        if ( obj != obj_end && *obj == '/' ) {
                            vt_i = 0;
                        } else {
                            if ( !parse_int( vt_i, obj, obj_end ) ) goto error;
                        }
                        dprint( "vt_i=" + std::to_string( vt_i ) );
                        vertex->vt_i = (vt_i >= 0) ? vt_i : ((int)hdr->texcoord_cnt + vt_i);

                        int vn_i;
                        if ( !expect_char( '/', obj, obj_end ) ) goto error;
                        if ( !parse_int( vn_i, obj, obj_end ) )  goto error ;
                        dprint( "vn_i=" + std::to_string( vn_i ) );
                        vertex->vn_i = (vn_i >= 0) ? vn_i : (hdr->norm_cnt + vn_i);

                    }
                    obj_assert( polygon->vtx_cnt != 0, ".obj f command has no vertices" );
                    if ( object != nullptr ) object->poly_cnt++;

                    // precompute surface normal and area (works for triangle only)
                    Vertex * pvertexes = &vertexes[polygon->vtx_i];
                    real3 p0( positions[pvertexes[0].v_i].c[0], positions[pvertexes[0].v_i].c[1], positions[pvertexes[0].v_i].c[2] );
                    real3 p1( positions[pvertexes[1].v_i].c[0], positions[pvertexes[1].v_i].c[1], positions[pvertexes[1].v_i].c[2] );
                    real3 p2( positions[pvertexes[2].v_i].c[0], positions[pvertexes[2].v_i].c[1], positions[pvertexes[2].v_i].c[2] );
                    polygon->normal = cross( sub( p1, p0 ), sub( p2, p0 ) );
                    real len = length( polygon->normal );
                    polygon->area = len / 2;
                    div( polygon->normal, len );
                    break;
                }

                case CMD_MTLLIB:
                    if ( !parse_name( mtllib, obj, obj_end ) ) goto error;
                    if ( !mtllib_load( dir_path, mtllib ) ) goto error;
                    break;
                    
                case CMD_USEMTL:
                    obj_assert( mtllib != nullptr, "no mtllib defined for object " + std::string( obj_name ) );
                    if ( !parse_name( name, obj, obj_end ) ) goto error;
                    mtl_name = std::string( name );
                    obj_assert( name_to_mtl.find( mtl_name ) != name_to_mtl.end(), "unknown material: " + std::string( mtl_name ) );
                    material = name_to_mtl[mtl_name];
                    mtl_i = material - materials;
                    break;

                default:
                    break;
            }
        }

    error:
        error_msg += " (at line " + std::to_string( line_num ) + " of " + obj_file + ")";
        assert( 0 );
    }

    ~Model() 
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

    bool write( std::string file_path, bool is_compressed=true ) 
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
            for( uint _byte_cnt = byte_cnt; _byte_cnt != 0;  ) \
            { \
                uint _this_byte_cnt = 1024*1024*1024; \
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

    bool write_uncompressed( std::string file_path )
    {
        int fd = open( file_path.c_str(), O_WRONLY );
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

        close( fd );
        return true;
    }

    Model( std::string file_path, bool is_compressed=true )
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
        // Write out header than individual arrays.
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
                for( uint _byte_cnt = (cnt)*sizeof(type); _byte_cnt != 0;  ) \
                { \
                    uint _this_byte_cnt = 1024*1024*1024; \
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
        _read( texels,      char,     hdr->texel_cnt );
        _read( strings,     char,     hdr->char_cnt );
        _read( bvh_nodes,   BVH_Node, hdr->bvh_node_cnt );

        gzclose( fd );

        is_good = true;
    }

    bool read_uncompressed( std::string file_path )
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
        _uread( texels,      char,     hdr->texel_cnt );
        _uread( strings,     char,     hdr->char_cnt );
        _uread( bvh_nodes,   BVH_Node, hdr->bvh_node_cnt );

        is_good = true;

        return true;
    }

private:
    typedef enum 
    {
        CMD_O,
        CMD_G,
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
        CMD_MAP_BUMP
    } mtl_cmd_t;
            
    char *              obj;
    char *              obj_start;
    char *              obj_end;
    uint                line_num;

    bool mtllib_load( std::string dir_path, std::string mtl_file )
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
        char * mtl_name = nullptr;
        char * tex_name = nullptr;

        Material * material = nullptr;
        Texture *  texture  = nullptr;
        uint       tex_i;

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
                    if ( !parse_name( mtl_name, mtl, mtl_end ) ) return false;
                    if ( name_to_mtl.find( mtl_name ) != name_to_mtl.end() ) {
                        material = name_to_mtl[ mtl_name ];         
                        dprint( "  found " + std::string( mtl_name ) );
                    } else {
                        perhaps_realloc<Material>( materials, hdr->mtl_cnt, max->mtl_cnt, 1 );
                        material = &materials[ hdr->mtl_cnt++ ];
                        memset( material, 0, sizeof( Material ) );
                        material->name_i = mtl_name - strings;
                        material->map_Ka_i = uint(-1);
                        material->map_Kd_i = uint(-1);
                        material->map_Ke_i = uint(-1);
                        material->map_Ks_i = uint(-1);
                        material->map_Bump_i = uint(-1);
                        name_to_mtl[ mtl_name ] = material;
                        dprint( "  added " + std::string( mtl_name ) );
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
                case CMD_MAP_BUMP:
                    rtn_assert( material != nullptr, "no material defined" );
                    if ( !parse_name( tex_name, mtl, mtl_end ) ) return false;
                    if ( name_to_tex.find( tex_name ) != name_to_tex.end() ) {
                        // already loaded it
                        //
                        texture = name_to_tex[ tex_name ];
                    } else {
                        // need to load it
                        //
                        if ( !load_texture( dir_path, tex_name, texture ) ) return false;
                    }

                    tex_i = texture - textures;

                    switch( cmd ) 
                    {
                        case CMD_MAP_KA:        material->map_Ka_i   = tex_i; break;
                        case CMD_MAP_KD:        material->map_Kd_i   = tex_i; break;
                        case CMD_MAP_KE:        material->map_Ke_i   = tex_i; break;
                        case CMD_MAP_KS:        material->map_Ks_i   = tex_i; break;
                        case CMD_MAP_BUMP:      material->map_Bump_i = tex_i; break;
                        default:                                              break; // should not happen
                    }

                    break;

                default:
                    rtn_assert( 0, "unknown .mtl command" );
            }

            rtn_assert( eol( mtl, mtl_end ), "not at eol in .mtl file" );
        }

        return true;
    }

    bool load_texture( std::string dir_path, char *& tex_name, Texture *& texture )
    {
        // allocate Texture structure
        //
        perhaps_realloc<Texture>( textures, hdr->tex_cnt, max->tex_cnt, 1 );
        texture = &textures[ hdr->tex_cnt++ ];
        memset( texture, 0, sizeof( Texture ) );
        texture->name_i = tex_name - strings;
        name_to_tex[ tex_name ] = texture;

        // convert file extension to .bmp (assume 3 characters in extension for now)
        // also change any '\' to '/'
        //
        char * bmp_name = strdup( tex_name );
        size_t len = strlen( bmp_name );
        for( size_t i = 0; i < len; i++ ) 
        {
            if ( bmp_name[i] == '\\' ) bmp_name[i] = '/';
        }
        char * ext = bmp_name + len - 4;
        rtn_assert( *ext == '.', "texture file does not have 3-character extension" );
        ext++;
        *ext = 'b'; ext++;
        *ext = 'm'; ext++;
        *ext = 'p'; ext++;

        // read .bmp file (the user is responsible for running mogrify in the textures directory)
        //
        char * data;
        char * data_end;
        if ( !open_and_read( dir_path, std::string( bmp_name ), data, data_end ) ) return false;

        // Extract image height and width from header.
        // Actual byte width is padded to 4-byte boundary.
        // Allocate RGB texels.
        //
        uint width  = *(int*)&data[18];
        uint height = *(int*)&data[22];
        texture->width  = width;
        texture->height = height;

        uint byte_width = width * 3;
        uint byte_cnt   = byte_width * height;

        char * bgr = data+54;  // first BGR texel
        texture->texel_i = hdr->texel_cnt;

        perhaps_realloc<char>( texels, hdr->texel_cnt, max->texel_cnt, byte_cnt );
        hdr->texel_cnt += byte_cnt;

        // Copy texels and convert BGR to RGB byte order.
        //
        char * rgb = texels + texture->texel_i;
        for( uint i = 0; i < height; i++ )
        {
            for( uint j = 0; j < byte_width; j += 3, bgr += 3, rgb += 3 )
            {
                if ( bgr < data_end ) {
                    rgb[0] = bgr[2];
                    rgb[1] = bgr[1];
                    rgb[2] = bgr[0];
                } else {
                    rgb[0] = 0;
                    rgb[1] = 0;
                    rgb[2] = 0;
                }
            }    
        }
        delete data;
        bgr = nullptr;

        // if requested, generate mipmap levels down to 1x1
        //
        while( hdr->mipmap_filter != MIPMAP_FILTER::NONE && !(width == 1 && height == 1) )
        {
            uint to_width  = width  >> 1;
            uint to_height = height >> 1;
            if ( to_width  == 0 ) to_width  = 1;
            if ( to_height == 0 ) to_height = 1;

            uint   to_byte_width = to_width * 3;
            uint   to_byte_cnt   = to_byte_width * to_height;
            perhaps_realloc<char>( texels, hdr->texel_cnt, max->texel_cnt, to_byte_cnt );
            hdr->texel_cnt += to_byte_cnt;

            rgb = texels + texture->texel_i;
            char * to_rgb        = rgb + byte_cnt;
            char * to_rgb_saved  = to_rgb;

            for( uint i = 0; i < to_height; i++ )
            {
                for( uint j = 0; j < to_byte_width; j += 3, to_rgb += 3 )
                {
                    uint sum[3] = { 0, 0, 0 };
                    uint cnt = 0;
                    for( uint fi = 2*i; fi < (2*i+1) && fi < height; fi++ )
                    {
                        char * frgb = rgb + fi*byte_width + 2*j;
                        for( uint fj = 2*j; fj < (2*j+1) && fj < byte_width; fj += 3, frgb += 3 )
                        {
                            sum[0] += frgb[0];
                            sum[1] += frgb[1];
                            sum[2] += frgb[2];
                            cnt++;
                        }
                    }

                    to_rgb[0] = (sum[0] + 0.5) / cnt;
                    to_rgb[1] = (sum[1] + 0.5) / cnt;
                    to_rgb[2] = (sum[2] + 0.5) / cnt;
                }
            }

            rgb        = to_rgb_saved;
            height     = to_height;
            width      = to_width;
            byte_width = to_byte_width;
        }

        return true;
    }

    bool open_and_read( std::string dir_path, std::string file_name, char *& start, char *& end )
    {
        std::string file_path = (dir_path != "") ? (dir_path + "/" + file_name) : file_name;
        const char * fname = file_path.c_str();
        int fd = open( fname, O_RDONLY );
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
            if ( read( fd, addr, _this_size ) <= 0 ) {
                close( fd );
                rtn_assert( 0, "could not read() file " + std::string(fname) + " - read error: " + std::string( strerror( errno ) ) );
            }
            size -= _this_size;
            addr += _this_size;
        }
        close( fd );
        return true;
    }

    inline void skip_whitespace( char *& xxx, char *& xxx_end )
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

    inline bool eol( char *& xxx, char *& xxx_end )
    {
        // skip some whitespace
        //
        while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;

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

    inline bool expect_char( char ch, char *& xxx, char* xxx_end )
    {
        rtn_assert( xxx != xxx_end, "premature end of file" );
        rtn_assert( *xxx == ch, "unexpected character: " + std::string( 1, *xxx ) );
        xxx++;
        return true;
    }

    inline bool expect_cmd( const char * s, char *& xxx, char *& xxx_end )
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

    inline bool parse_name( char *& name, char *& xxx, char *& xxx_end )
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

        rtn_assert( 0, "could not parse name" );
    }

    inline bool parse_obj_cmd( obj_cmd_t& cmd )
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
                rtn_assert( 0, "bad .obj command character: " + std::string( 1, ch ) );
        }
    }

    inline bool parse_mtl_cmd( mtl_cmd_t& cmd, char *& mtl, char *& mtl_end )
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

                rtn_assert( 0, "truncated .mtl m command" );

            default:
                rtn_assert( 0, "bad .mtl command character" );
        }
    }

    inline bool parse_real3( real3& r3, char *& xxx, char *& xxx_end )
    {
        return parse_real( r3.c[0], xxx, xxx_end ) && 
               parse_real( r3.c[1], xxx, xxx_end ) && 
               parse_real( r3.c[2], xxx, xxx_end );
    }

    inline bool parse_real2( real2& r2, char *& xxx, char *& xxx_end )
    {
        return parse_real( r2.c[0], xxx, xxx_end ) && 
               parse_real( r2.c[1], xxx, xxx_end );
    }

    inline bool parse_real( real& r, char *& xxx, char *& xxx_end )
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
        rtn_assert( vld, "unable to parse real in " + std::string( ((xxx == obj) ? ".obj" : ".mtl") ) + " file" );
        return vld;
    }

    inline bool parse_int( int& i, char *& xxx, char *& xxx_end )
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

    inline real dot( const real3 &v1, const real3 &v2 ) 
    {
        return v1.c[0] *v2.c[0] + v1.c[1] *v2.c[1]  + v1.c[2] *v2.c[2];
    }

    inline real3 cross( const real3 &v1, const real3 &v2 ) 
    {
        return real3( (v1.c[1]*v2.c[2] - v1.c[2]*v2.c[1]),
                      (-(v1.c[0]*v2.c[2] - v1.c[2]*v2.c[0])),
                      (v1.c[0]*v2.c[1] - v1.c[1]*v2.c[0]) );
    }

    inline real length( const real3& v ) const 
    { 
        return std::sqrt( v.c[0]*v.c[0] + v.c[1]*v.c[1] + v.c[2]*v.c[2] ); 
    }

    inline real3 add( const real3& v1, const real3& v2 )
    {
        real3 r;
        r.c[0] = v1.c[0] + v2.c[0];
        r.c[1] = v1.c[1] + v2.c[1];
        r.c[2] = v1.c[2] + v2.c[2];
        return r;
    }

    inline real3 sub( const real3& v1, const real3& v2 )
    {
        real3 r;
        r.c[0] = v1.c[0] - v2.c[0];
        r.c[1] = v1.c[1] - v2.c[1];
        r.c[2] = v1.c[2] - v2.c[2];
        return r;
    }

    inline void div( real3& v, real denom ) 
    {
        v.c[0] /= denom;
        v.c[1] /= denom;
        v.c[2] /= denom;
    }

    void bvh_build( BVH_TREE bvh_tree )
    {
        (void)bvh_tree;
        hdr->bvh_root_i = bvh_node( 0, hdr->poly_cnt, 1 );
    }

    void poly_bounding_box( uint poly_i, AABB& box )
    {
        const Polygon& poly = polygons[poly_i];
        for( uint i = 0; i < poly.vtx_cnt; i++ )
        {
            const Vertex& vtx = vertexes[poly.vtx_i+i];
            const real3&  pos = positions[vtx.v_i];
            if ( i == 0 ) {
                box.min = pos;
                box.max = pos;
            } else {
                for( uint j = 0; j < 3; j++ )
                {
                    if ( pos.c[j] < box.min.c[j] ) box.min.c[j] = pos.c[j];
                    if ( pos.c[j] > box.max.c[j] ) box.max.c[j] = pos.c[j];
                }
            }
        }
    }

    inline uint bvh_qsplit( uint poly_i, uint n, real pivot, uint axis )
    {
       uint m = poly_i;

       for( uint i = poly_i; i < (poly_i+n); i++ )
       {
           AABB box;
           poly_bounding_box( i, box );
           real centroid = (box.min.c[axis] + box.max.c[axis]) * 0.5;
           if ( centroid < pivot ) {
               Polygon temp = polygons[i];
               polygons[i]  = polygons[m];
               polygons[m]  = temp;
               m++;
           }
        }

        if ( m == poly_i || m == (poly_i+n) ) m = poly_i + n/2;
        return m;
    }

    uint bvh_node( uint poly_i, uint n, uint axis ) 
    {
        perhaps_realloc( bvh_nodes, hdr->bvh_node_cnt, max->bvh_node_cnt, 1 );
        uint bvh_i = hdr->bvh_node_cnt++;
        BVH_Node * node = &bvh_nodes[bvh_i];

        if ( n == 1 || n == 2 ) {
            node->left_is_leaf = true;
            node->right_is_leaf = true;
            node->left_i  = poly_i;
            node->right_i = poly_i + n-1;
            poly_bounding_box( poly_i, node->box );
            if ( n == 2 ) {
                AABB new_box;
                poly_bounding_box( poly_i+1, new_box );
                node->box.expand( new_box );
            }

        } else {
            node->left_is_leaf = false;
            node->right_is_leaf = false;
            real pivot = (node->box.min.c[axis] + node->box.max.c[axis]) * 0.5;
            uint m = bvh_qsplit( poly_i, n, pivot, axis );
            uint nm = m - poly_i + 1;
            node->left_i  = bvh_node( poly_i, nm,   (axis + 1) % 3 );
            node->right_i = bvh_node(      m, n-nm, (axis + 1) % 3 );

            node->box = bvh_nodes[node->left_i].box;
            node->box.expand( bvh_nodes[node->right_i].box );
        }

        return bvh_i;
    }
};

#endif
