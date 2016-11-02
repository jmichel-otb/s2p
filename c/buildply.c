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

void write_ply_header(FILE* f, uint64_t npoints, int zone,
        bool hem, bool colors)
{
    fprintf(f, "ply\n");
    fprintf(f, "format binary_little_endian 1.0\n");
    fprintf(f, "comment created by S2P\n");
    if (zone >= 0)
        fprintf(f, "comment projection: UTM %i%s\n", zone, (hem ? "N" : "S"));
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
    fprintf(stderr, "\t usage: %s out.ply ecef_coord_crop.tif [roi_color_ref_crop.tif]", s);
}

int main(int c, char *v[])
{

    if (!( (c==3) || (c==4) ) ) {
        help(*v);
        return 1;
    }
    
    bool there_is_color=false;
    if (c==4)
        there_is_color = true;
    
    // parse the remaining arguments
    char *fname_ply = v[1];
    char *fname_ecef = v[2];
    char *fname_colors;
    if (there_is_color)
        fname_colors = v[3];

    // read input images
    int w, h, pd;
    double *ecef = iio_read_image_double_vec(fname_ecef, &w, &h, &pd);
    int wc, hc, pdc;
    double *clr;
    if (there_is_color)
         clr = iio_read_image_double_vec(fname_colors, &wc, &hc, &pdc);
         
    if ( (w != wc) || (h != hc) ) 
    { 
        printf("color and ecef image size mismatch\n");
        return 1;
    }
    if (pd != 3) 
    {
        printf("ecef image must have 3 bands\n");
        return 1;
    }
        
    int zone=-1;
    bool hem;

    // count number of valid pixels
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
    write_ply_header(ply_file, npoints, zone, hem, 
                    there_is_color);

    // fill ply file
    size_t point_size = 3*sizeof(double);
    int dim = 3;
    if (there_is_color)
    {
        point_size += 3*sizeof(uint8_t);
        dim+=3;
    }
    
    char *buf = (char *) malloc(point_size);
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
            
            // write to memory
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
        }
    }

    free(buf);
    if (there_is_color)
        free(clr);
    free(utm_coord);
    free(ecef);
    fclose(ply_file);
    
    return 0;
}
