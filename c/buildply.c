#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "coordconvert.h"

void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);

static void get_utm_coord(double *out, double lat, double lon, double alt, int zone)
{
    utm_alt_zone(out, lat, lon, zone);
    out[2] = alt;
}

void write_ply_header(FILE* f, uint64_t npoints, char * comments, bool colors)
{
    fprintf(f, "ply\n");
    fprintf(f, "format binary_little_endian 1.0\n");
    fprintf(f, "comment created by S2P\n");
    fprintf(f, "%s",comments);
    /*if (zone >= 0)
        fprintf(f, "comment projection: UTM %i%s\n", zone, (hem ? "N" : "S"));*/
    fprintf(f, "element vertex %" PRIu64 "\n", npoints);
    fprintf(f, "property double x\n");
    fprintf(f, "property double y\n");
    fprintf(f, "property double z\n");
    if (colors) {
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
    }
    fprintf(f, "end_header\n");
}


void help(char *s)
{
    fprintf(stderr, "\t usage: %s root_out_dir /*out.ply ecef_coord_crop.tif [roi_color_ref_crop.tif]*/", s);
}

int main(int c, char *v[])
{

    if ( c!= 2 ) {
        help(*v);
        return 1;
    }
    
    // parse arguments
    char *root_out_dir=v[1];
    char fname_ply[1000];
    sprintf(fname_ply,"%s/cloud.ply",root_out_dir);
    char fname_ecef_ply[1000];
    sprintf(fname_ecef_ply,"%s/cloud_ecef.ply",root_out_dir);
    char fname_ecef[1000];
    sprintf(fname_ecef,"%s/ecef_coord_crop.tif",root_out_dir);
    char fname_colors[1000];
    sprintf(fname_colors,"%s/roi_color_ref_crop.tif",root_out_dir);
    
    bool there_is_color=false;
    FILE *file;
    if (file = fopen(fname_colors, "r"))
    {
        fclose(file);
        there_is_color = true;
    }

    // read input images
    // * ecef coord
    int w, h, pd;
    double *ecef = iio_read_image_double_vec(fname_ecef, &w, &h, &pd);
    // * colors
    int wc, hc, pdc;
    double *clr;
    if (there_is_color)
         clr = iio_read_image_double_vec(fname_colors, &wc, &hc, &pdc);
         
    if ( there_is_color && ( (w != wc) || (h != hc) ) )
    { 
        printf("color and ecef image size mismatch\n");
        return 1;
    }
    if (pd != 3) 
    {
        printf("ecef image must have 3 bands\n");
        return 1;
    }
        

    // count number of valid pixels
    int zone=-1;
    bool hem;
    uint64_t npoints = 0;
    printf("counting valid points...\n");
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) 
    {
        //uint64_t pix = (uint64_t) row * w + col;
        bool ok=true;
        for(int t=0; t<3; t++)
            if (isnan(ecef[w*3*y+3*x+t]))
                ok=false;
        
        if (ok) 
        {
            npoints++;

            // UTM Zone will be the zone of first 'not NaN' point
            if (zone < 0) 
            {
                double lgt,lat,alt;
                ECEF_to_lgt_lat_alt(ecef[w*3*y+3*x], 
                                    ecef[w*3*y+3*x+1], 
                                    ecef[w*3*y+3*x+2],
                                    &lgt,&lat,&alt);
                utm_zone(&zone, &hem, lat, lgt);
            }
        }
    }
    printf("found %" PRIu64 " valid points\n", npoints);

    // print header for ply file
    FILE *ply_file = fopen(fname_ply, "wb");
    char comments_UTM[1000];
    sprintf(comments_UTM,"comment projection: UTM %i%s\n", zone, (hem ? "N" : "S"));
    write_ply_header(ply_file, npoints, comments_UTM, there_is_color);
    
    FILE *ecef_ply_file = fopen(fname_ecef_ply, "wb");
    char comments_ECEF[1000];
    sprintf(comments_ECEF,"comment projection: ECEF\n");
    write_ply_header(ecef_ply_file, npoints, comments_ECEF, there_is_color);

    // fill ply file
    size_t point_size = 3*sizeof(double);
    int dim = 3;
    if (there_is_color)
    {
        point_size += 3*sizeof(uint8_t);
        dim+=3;
    }
    
    char *buf = (char *) malloc(point_size);
    char *buf_ecef = (char *) malloc(point_size);
    double *utm_coord = (double *) malloc(dim*sizeof(double));
    
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) 
    {
        bool ok=true;
        for(int t=0; t<3; t++)
            if (isnan(ecef[w*3*y+3*x+t]))
                ok=false;
        
        if (ok) 
        {
            // positions
            double lgt,lat,alt;
            ECEF_to_lgt_lat_alt(ecef[w*3*y+3*x], 
                                ecef[w*3*y+3*x+1], 
                                ecef[w*3*y+3*x+2],
                                &lgt,&lat,&alt);
            get_utm_coord(utm_coord, lat, lgt, alt, zone);
            
            // colors
            if (there_is_color)
                for(int t=0; t<3; t++)
                    utm_coord[t+3] = clr[w*3*y+3*x+t];
            
            // write to memory UTM coord
            double *ptr_double = (double *) buf;
            ptr_double[0] = utm_coord[0];
            ptr_double[1] = utm_coord[1];
            ptr_double[2] = utm_coord[2];
            
            if (there_is_color)
            {
                char *ptr_char = buf + 3*sizeof(double);
                ptr_char[0] = utm_coord[3];
                ptr_char[1] = utm_coord[4];
                ptr_char[2] = utm_coord[5];
            }
            fwrite( buf, point_size, 1, ply_file);
            
            // write to memory ECEF coord
            ptr_double = (double *) buf_ecef;
            ptr_double[0] = ecef[w*3*y+3*x];
            ptr_double[1] = ecef[w*3*y+3*x+1];
            ptr_double[2] = ecef[w*3*y+3*x+2];
            
            if (there_is_color)
            {
                char *ptr_char = buf_ecef + 3*sizeof(double);
                ptr_char[0] = utm_coord[3];
                ptr_char[1] = utm_coord[4];
                ptr_char[2] = utm_coord[5];
            }
            fwrite( buf_ecef, point_size, 1, ecef_ply_file);
        }
    }

    free(buf);
    free(buf_ecef);
    if (there_is_color)
        free(clr);
    free(utm_coord);
    free(ecef);
    fclose(ply_file);
    fclose(ecef_ply_file);
    
    return 0;
}
