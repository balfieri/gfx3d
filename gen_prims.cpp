// gen_prims.cpp - generates prims from .obj format models
//
static inline float uniform()  { return 0.0; }    // TODO: temporary

#include "Model.h"
#include <iostream>
#include <stdlib.h>

static inline void print( std::string msg )
{
    std::cout << msg << "\n";
}

static inline void die( std::string msg )             
{ 
    print( "ERROR: " + msg );
    exit( 1 );
}

int main( int argc, const char * argv[] )
{
    print( "" );
    std::string file_name;
    if ( argc == 2 ) {
        file_name = std::string( argv[1] );
    } else {
        die( "usage: gen_prims <top_file>" );
    }
    std::cout << "Reading in " << file_name << "...\n";
    Model * model = new Model( file_name, Model::MIPMAP_FILTER::BOX, Model::BVH_TREE::BINARY );
    if ( !model->is_good ) die( model->error_msg );

    print( "char_cnt:       " + std::to_string( model->hdr->char_cnt ) );
    print( "obj_cnt:        " + std::to_string( model->hdr->obj_cnt ) );
    print( "poly_cnt:       " + std::to_string( model->hdr->poly_cnt ) );
    print( "vtx_cnt:        " + std::to_string( model->hdr->vtx_cnt ) );
    print( "pos_cnt:        " + std::to_string( model->hdr->pos_cnt ) );
    print( "norm_cnt:       " + std::to_string( model->hdr->norm_cnt ) );
    print( "texcoord_cnt:   " + std::to_string( model->hdr->texcoord_cnt ) );
    print( "mtl_cnt:        " + std::to_string( model->hdr->mtl_cnt ) );
    print( "tex_cnt:        " + std::to_string( model->hdr->tex_cnt ) );
    print( "texel_cnt:      " + std::to_string( model->hdr->texel_cnt ) );
    print( "bvh_node_cnt:   " + std::to_string( model->hdr->bvh_node_cnt ) );
    print( "matrix_cnt:     " + std::to_string( model->hdr->matrix_cnt ) );
    print( "inst_cnt:       " + std::to_string( model->hdr->inst_cnt ) );
    print( "light_cnt:      " + std::to_string( model->hdr->light_cnt ) );
    print( "camera_cnt:     " + std::to_string( model->hdr->camera_cnt ) );
    print( "frame_cnt:      " + std::to_string( model->hdr->frame_cnt ) );
    print( "animation_cnt:  " + std::to_string( model->hdr->animation_cnt ) );
    float mb_cnt = float(model->hdr->byte_cnt) / float(1024 * 1024);
    print( "total_byte_cnt: " + std::to_string( model->hdr->byte_cnt ) + " (" + std::to_string( mb_cnt ) + " MB)" );

    // construct name of .model file 
    std::string dir_name;
    std::string base_name;
    std::string ext_name;
    Model::dissect_path( file_name, dir_name, base_name, ext_name );
    if ( dir_name == std::string( "" ) ) die( "top_file has no dir_name: " + file_name );
    std::string dir_name2;
    Model::dissect_path( dir_name, dir_name2, base_name, ext_name );
    if ( dir_name2 == std::string( "" ) ) die( "something really bad happened: dir_name=" + dir_name + 
                                               " dir_name2=" + dir_name2 + " base_name=" + base_name + " ext_name=" + ext_name );
    std::string name = dir_name + "/" + dir_name2 + ".model";

    Model * model2;
    if ( 0 ) {
        std::string gz_name = name;
        print( "Writing compressed " + gz_name + "..." );
        if ( !model->write( gz_name, true ) ) die( model->error_msg );

        print( "Reading compressed " + gz_name + "..." );
        model2 = new Model( gz_name, true );
        if ( !model2->is_good ) die( model2->error_msg );
        if ( model2->hdr->obj_cnt != model->hdr->obj_cnt ||
             model2->hdr->poly_cnt != model->hdr->poly_cnt ||
             model2->hdr->vtx_cnt != model->hdr->vtx_cnt ||
             model2->hdr->pos_cnt != model->hdr->pos_cnt ||
             model2->hdr->norm_cnt != model->hdr->norm_cnt ||
             model2->hdr->texcoord_cnt != model->hdr->texcoord_cnt ||
             model2->hdr->mtl_cnt != model->hdr->mtl_cnt ||
             model2->hdr->tex_cnt != model->hdr->tex_cnt ||
             model2->hdr->texel_cnt != model->hdr->texel_cnt ||
             model2->hdr->char_cnt != model->hdr->char_cnt ) {
            die( "read in model's hdr does not match previous" );
        }
        delete model2;
    }

    print( "Writing uncompressed " + name + "..." );
    if ( !model->write( name, false ) ) die( model->error_msg );

    print( "Reading uncompressed " + name + "..." );
    model2 = new Model( name, false );
    if ( !model2->is_good ) die( model2->error_msg );
    if ( model2->hdr->obj_cnt != model->hdr->obj_cnt ||
         model2->hdr->poly_cnt != model->hdr->poly_cnt ||
         model2->hdr->vtx_cnt != model->hdr->vtx_cnt ||
         model2->hdr->pos_cnt != model->hdr->pos_cnt ||
         model2->hdr->norm_cnt != model->hdr->norm_cnt ||
         model2->hdr->texcoord_cnt != model->hdr->texcoord_cnt ||
         model2->hdr->mtl_cnt != model->hdr->mtl_cnt ||
         model2->hdr->tex_cnt != model->hdr->tex_cnt ||
         model2->hdr->texel_cnt != model->hdr->texel_cnt ||
         model2->hdr->char_cnt != model->hdr->char_cnt ||
         model2->hdr->bvh_node_cnt != model->hdr->bvh_node_cnt ) {
        die( "read in model's hdr does not match previous" );
    }

    delete model;
    delete model2;

    print( "PASSED" );
    return 0;
}
