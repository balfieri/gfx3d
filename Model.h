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
//              mogrify *.{jpg,png} -format bmp
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
//     3) This code is set up to write out everything to one binary file, but I don't think we need to do that
//        yet given how fast this is.
//
#ifndef _Model_h
#define _Model_h

#include <cstdint>
#include <string>
#include <map>
#include <iostream>

#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

class Model
{
public:
    typedef uint32_t uint;
    typedef float    real;

    class real3
    {
    public:
        real            c[3];
    };

    class real2
    {
    public:
        real            c[2];
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
        uint        mtl_cnt;                // in materials array
        uint        tex_cnt;                // in textures array
        uint        texel_cnt;              // in texels array  (last in file)
        uint        char_cnt;               // in strings array
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
        uint        name_i;                 // index in strings array (null-terminated strings)
        uint        width;
        uint        height;
        uint        texel_i;                // index into texels array of first texel
    };

    // public fields
    //
    static const uint VERSION = 0xB0BA1f01; // current version is 1

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    // structs
    Header              hdr;
    Header              max;                // holds max lengths of currently allocated arrays 

    // arrays
    Object *            objects;
    Polygon *           polygons;
    Vertex *            vertexes;
    real3 *             positions;
    real3 *             normals;
    real2 *             texcoords;
    Material *          materials;
    Texture *           textures;
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
    #define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = msg; return false; }
    #define obj_assert( bool, msg ) if ( !(bool) ) { error_msg = msg; goto error;   }

    Model( std::string dir_path, std::string obj_file )
    {
        is_good = false;

        memset( &hdr, 0, sizeof( Header ) );
        hdr.version = VERSION;
        hdr.pos_cnt = 1;
        hdr.norm_cnt = 1;
        hdr.texcoord_cnt = 1;

        //------------------------------------------------------------
        // Initial lengths of arrays are large in virtual memory
        //------------------------------------------------------------
        max.obj_cnt     =  1000000;
        max.poly_cnt    = 20000000;
        max.mtl_cnt     =     1000;

        max.vtx_cnt     = 3*max.poly_cnt;
        max.pos_cnt     = max.poly_cnt;
        max.norm_cnt    = max.poly_cnt;
        max.texcoord_cnt= max.poly_cnt;
        max.tex_cnt     = max.mtl_cnt;
        max.texel_cnt   = max.mtl_cnt * 1000000;
        max.char_cnt    = (max.obj_cnt + max.mtl_cnt) * 128;

        //------------------------------------------------------------
        // Allocate arrays
        //------------------------------------------------------------
        objects         = new Object[   max.obj_cnt ];
        polygons        = new Polygon[  max.poly_cnt ];
        vertexes        = new Vertex[   max.vtx_cnt ];
        positions       = new real3[    max.pos_cnt ];
        normals         = new real3[    max.norm_cnt ];
        texcoords       = new real2[    max.texcoord_cnt ];
        materials       = new Material[ max.mtl_cnt ];
        textures        = new Texture[  max.tex_cnt ];
        texels          = new char[     max.texel_cnt ];
        strings         = new char[     max.char_cnt ];

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
        Material *  material;
        uint        mtl_i = uint(-1);
        Object *    object;
        Polygon *   polygon;
        Vertex *    vertex;

        for( ;; ) 
        {
            skip_whitespace( obj, obj_end );
            if ( obj == obj_end ) {
                // done, no errors
                hdr.byte_cnt = uint64_t( 1                ) * sizeof( hdr ) +
                               uint64_t( hdr.obj_cnt      ) * sizeof( objects[0] ) +
                               uint64_t( hdr.poly_cnt     ) * sizeof( polygons[0] ) +
                               uint64_t( hdr.vtx_cnt      ) * sizeof( vertexes[0] ) +
                               uint64_t( hdr.pos_cnt      ) * sizeof( positions[0] ) +
                               uint64_t( hdr.norm_cnt     ) * sizeof( normals[0] ) +
                               uint64_t( hdr.texcoord_cnt ) * sizeof( texcoords[0] ) +
                               uint64_t( hdr.mtl_cnt      ) * sizeof( materials[0] ) +
                               uint64_t( hdr.tex_cnt      ) * sizeof( textures[0] ) +
                               uint64_t( hdr.texel_cnt    ) * sizeof( texels[0] ) +
                               uint64_t( hdr.char_cnt     ) * sizeof( strings[0] );
                is_good = true;
                return;
            }

            obj_cmd_t cmd;
            if ( !parse_obj_cmd( cmd ) ) break;
            dprint( "obj_cmd=" + std::to_string( cmd ) );

            switch( cmd )
            {
                case CMD_O:
                    obj_assert( hdr.obj_cnt < max.obj_cnt, "objects[] is full" );
                    object = &objects[ hdr.obj_cnt++ ];

                    if ( !parse_name( obj_name, obj, obj_end ) ) goto error;
                    object->name_i = obj_name - strings;
                    object->poly_cnt = 0;
                    object->poly_i = hdr.poly_cnt;
                    break;
                    
                case CMD_G:
                    if ( !parse_name( name, obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_V:
                    obj_assert( hdr.pos_cnt < max.pos_cnt, "positions[] is full" );
                    if ( !parse_real3( positions[ hdr.pos_cnt++ ], obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_VN:
                    obj_assert( hdr.norm_cnt < max.norm_cnt, "normals[] is full" );
                    if ( !parse_real3( normals[ hdr.norm_cnt++ ], obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_VT:
                    obj_assert( hdr.texcoord_cnt < max.texcoord_cnt, "texcoords[] is full" );
                    if ( !parse_real2( texcoords[ hdr.texcoord_cnt++ ], obj, obj_end ) ) goto error;
                    break;
                    
                case CMD_F:
                    obj_assert( hdr.poly_cnt < max.poly_cnt, "polygons[] is full" );
                    polygon = &polygons[ hdr.poly_cnt++ ];
                    polygon->mtl_i = mtl_i;
                    polygon->vtx_cnt = 0;
                    polygon->vtx_i = hdr.vtx_cnt;
                    while( !eol( obj, obj_end ) ) 
                    {
                        polygon->vtx_cnt++;
                        obj_assert( hdr.vtx_cnt < max.vtx_cnt, "vertexes[] is full" );
                        vertex = &vertexes[ hdr.vtx_cnt++ ];

                        int v_i;
                        if ( !parse_int( v_i, obj, obj_end ) )   goto error;
                        dprint( "v_i=" + std::to_string( v_i ) );
                        vertex->v_i = (v_i >= 0)  ? v_i : (hdr.pos_cnt + v_i);

                        if ( !expect_char( '/', obj, obj_end ) ) goto error;

                        int vt_i;
                        if ( obj != obj_end && *obj == '/' ) {
                            vt_i = 0;
                        } else {
                            if ( !parse_int( vt_i, obj, obj_end ) ) goto error;
                        }
                        dprint( "vt_i=" + std::to_string( vt_i ) );
                        vertex->vt_i = (vt_i >= 0) ? vt_i : ((int)hdr.texcoord_cnt + vt_i);

                        int vn_i;
                        if ( !expect_char( '/', obj, obj_end ) ) goto error;
                        if ( !parse_int( vn_i, obj, obj_end ) )  goto error ;
                        dprint( "vn_i=" + std::to_string( vn_i ) );
                        vertex->vn_i = (vn_i >= 0) ? vn_i : (hdr.norm_cnt + vn_i);

                    }
                    object->poly_cnt++;
                    obj_assert( polygon->vtx_cnt != 0, ".obj f command has no vertices" );
                    break;

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

            if ( !eol( obj, obj_end ) ) {
                error_msg += " in .mtl file";
                goto error;
            }
        }

    error:
        error_msg += " (at line " + std::to_string( line_num ) + " of " + obj_file + ")";
    }

    ~Model() {}

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
                        rtn_assert( hdr.mtl_cnt < max.mtl_cnt, "materials[] is full" );
                        material = &materials[ hdr.mtl_cnt++ ];
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
                    return false;
            }

            if ( !eol( mtl, mtl_end ) ) {
                error_msg += " in .mtl file";
                return false;
            }
        }

        return true;
    }

    bool load_texture( std::string dir_path, char *& tex_name, Texture *& texture )
    {
        // allocate Texture structure
        //
        rtn_assert( hdr.tex_cnt < max.tex_cnt, "textures[] is full" );
        texture = &textures[ hdr.tex_cnt++ ];
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
        uint actual_byte_width = (byte_width + 3) & ~0x3;
        uint actual_byte_cnt   = actual_byte_width * height;
        uint pad_byte_width    = actual_byte_width - byte_width;

        char * bgr = data+54;  // first BGR texel
        char * rgb = texels + hdr.texel_cnt;
        hdr.texel_cnt += actual_byte_cnt;
        rtn_assert( hdr.texel_cnt <= max.texel_cnt, "ran out of texel space" );

        // Copy texels and convert BGR to RGB byte order.
        //
        for( uint i = 0; i < height; i++ )
        {
            for( uint j = 0; j < byte_width; j += 3, bgr += 3, rgb += 3 )
            {
                rgb[0] = bgr[2];
                rgb[1] = bgr[1];
                rgb[2] = bgr[0];
            }    
            for( uint j = 0; j < pad_byte_width; j++, bgr++, rgb++ )
            {
                *rgb = *bgr;  // should be zeros, but copy anyway
            }
        }

        return true;
    }

    bool open_and_read( std::string dir_path, std::string file_name, char *& start, char *& end )
    {
        std::string file_path = dir_path + "/" + file_name;
        const char * fname = file_path.c_str();
        int fd = open( fname, O_RDONLY );
        rtn_assert( fd >= 0, "could not open file " + file_path + " - open() error: " + strerror( errno ) );

        struct stat file_stat;
        int status = fstat( fd, &file_stat );
        if ( status < 0 ) {
            close( fd );
            error_msg = "could not stat file " + std::string(fname) + " - stat() error: " + strerror( errno );
            return false;
        }
        size_t size = file_stat.st_size;

        // this large read should behave like an mmap() inside the o/s kernel and be as fast
        start = new char[ size ];
        if ( start == nullptr ) {
            close( fd );
            error_msg = "could not read file " + std::string(fname) + " - malloc() error: " + strerror( errno );
            return false;
        }
        if ( read( fd, start, size ) <= 0 ) {
            close( fd );
            error_msg = "could not read() file " + std::string(fname) + " - read error: " + std::string( strerror( errno ) );
            return false;
        }
        close( fd );
        end = start + size;
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
            return true;
        } else {
            error_msg = "end-of-line is not next, found character: " + std::string( 1, *xxx );
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
        name = &strings[hdr.char_cnt];

        while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

        uint len = 0;
        while( xxx != xxx_end )
        {
            char ch = *xxx;
            if ( ch == '\n' || ch == '\r' ) break;

            name[len++] = ch;
            vld = true;
            xxx++;
        }

        if ( vld ) {
            name[len] = '\0';
            hdr.char_cnt += len+1;
            char * ptr;
            for( ptr = &name[len-1]; ptr != name; ptr-- )
            {
                // skip trailing spaces
                if ( *ptr != ' ' && *ptr != '\t' ) break;
                *ptr = '\0';
            }

            return *ptr != '\0';
        }

        error_msg = "could not parse name";
        return false;
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

                error_msg = "bad .obj V command";
                return false;

            case 'm':
                cmd = CMD_MTLLIB;
                return expect_cmd( "tllib", obj, obj_end ); 
                
            case 'u':
                cmd = CMD_USEMTL;
                return expect_cmd( "semtl", obj, obj_end ); 
                
            default:
                error_msg = "bad .obj command character: " + std::string( 1, ch );
                return false;
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

                error_msg = "bad .mtl N command";
                return false;

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

                error_msg = "truncated .mtl K command";
                return false;

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

                error_msg = "truncated .mtl T command";
                return false;

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
                        error_msg = "bad .mtl map_k command"; 
                        return false;
                    }

                    mtl++;
                    if ( !expect_char( ' ', mtl, mtl_end ) ) {
                        error_msg = "bad .mtl map_k command (no space after Ka/Kd/Ke/Ks)";
                        return false;
                    }
                    return true;

                } else if ( ch == 'b' || ch == 'B' ) {
                    if ( !expect_char( 'u', mtl, mtl_end ) || 
                         !expect_char( 'm', mtl, mtl_end ) || 
                         !expect_char( 'p', mtl, mtl_end ) ||
                         !expect_char( ' ', mtl, mtl_end ) ) {
                        error_msg = "unexpected .mtl map_b command";
                        return false;
                    }

                    cmd = CMD_MAP_BUMP;
                    return true;
                }

                error_msg = "truncated .mtl m command";
                return false;

            default:
                error_msg = "bad .mtl command character";
                return false;
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
        if ( !vld ) error_msg = "unable to parse real in " + std::string( ((xxx == obj) ? ".obj" : ".mtl") ) + " file";
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
};

#endif
