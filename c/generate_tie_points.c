#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"
#include "coordconvert.h"


int main_generate_tie_points(int c, char *v[])
{
    if (c != 7) {
        fprintf(stderr, "usage:\n\t"
                "%s global_out_dir, col, row, tw, th, z "
              // 0         1         2    3   4   5   6      
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }
    

    // read input data
    char *global_out_dir=v[1];
    int X=atoi(v[2]);
    int Y=atoi(v[3]);
    int width=atoi(v[4]);
    int height=atoi(v[5]);
    double z = atof(v[6]);
    
    // build tile_dir path
    char tile_dir[1000];
    sprintf(tile_dir,"%s/tile_%d_%d_row_%d/col_%d",global_out_dir,width,height,Y,X);
    
    // take into account zoom
    width = (int) (width/z)+1;
    height = (int) (height/z)+1;
    
    // ##################################
    // build outputs paths and alloc mem
    // ##################################
    
    // * ecef
    char fout_ecef[1000];
    sprintf(fout_ecef,"%s/ecef_coord.tif",tile_dir);
    //float *ecef = (float *) calloc(width*height, sizeof(float));
    int wid,hei,pd;
    float *ecef = NULL;
    ecef = iio_read_image_float_vec(fout_ecef, &wid, &hei,&pd);
    
    if (!ecef)
        return -1;
    
    if ( (wid != width) || (hei != height) || (pd != 3) )
        return -2;

    
    // Yet another interesting thing :
    // output coord of the center of the
    // tile to eventually refine rpc coefs
    int mid_height = height / 2;
    int mid_width = width / 2;
    double dist_min=1e9;
    int best_x,best_y;
    double best_mid_height, best_mid_width;
    double best_lgt,best_lat,best_alt;
    char fname_tie_points[1000];
    bool refine_rpc=false;
    
    
    // ################################
    // time to get some tie points
    // ################################
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
        {
            // position inside ecef map
            if ( !isnan(ecef[width*3*y+3*x+0]) || !isnan(ecef[width*3*y+3*x+1]) || !isnan(ecef[width*3*y+3*x+2]) )
            {
                double global_x = x*z+X;
                double global_y = y*z+Y;
            
                double lgt,lat,alt;
                
                ECEF_to_lgt_lat_alt(ecef[width*3*y+3*x+0], 
                                ecef[width*3*y+3*x+1], 
                                ecef[width*3*y+3*x+2],
                                    &lgt,&lat,&alt);

                
                // Remember coord of the center of the tile
                // to eventually refine rpc
                double dist =  pow(x - mid_width,2.0);
                dist += pow(y - mid_height,2.0);
                dist = sqrt(dist);
                if ( (dist<dist_min) && (!isnan(lgt)) && (!isnan(lat)) && (!isnan(alt)) )
                {
                  dist_min=dist;
                  best_x = x;
                  best_y = y;
                  best_mid_height = global_y;
                  best_mid_width = global_x;
                  best_lgt = lgt;
                  best_lat = lat;
                  best_alt = alt;
                  refine_rpc=true;
                }
            }
        }

    // clean mem
    free(ecef);
    
    // rpc refining
    FILE * ftie_points = NULL;
    if ( refine_rpc )
    {
        sprintf(fname_tie_points,"%s/tie_points_bis.txt",tile_dir);
        ftie_points = fopen(fname_tie_points,"w");
        if (ftie_points)
        {
            fprintf(ftie_points, "%.2f %.2f %.6f %.6f %.3f\n",
                best_mid_width,best_mid_height,best_lgt,best_lat,best_alt);
            fclose(ftie_points);
        }
    }
    
    return 0;
}

int main(int c, char *v[])
{
    return main_generate_tie_points(c, v);
}
