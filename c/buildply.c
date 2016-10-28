#include <inttypes.h>
#include <stdio.h>
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
        bool hem, bool colors, bool normals)
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

int main(int c, char *v[])
{

    // parse the remaining arguments
    char *fname_ply = v[1];
    char *fname_ecef = v[2];

    // read input images
    int w, h, pd;
    double *ecef = iio_read_image_double_vec(fname_ecef, &w, &h, &pd);
    //printf("%d %d %d\n",w,h,pd);
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
                    false, false);

    // fill ply file
    size_t point_size = 3;
    double *utm_coord = (double *) malloc(point_size*sizeof(double));
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) 
    {
        bool ok=true;
        for(int t=0; t<3; t++)
            if (isnan(ecef[w*3*y+3*x+t]))
                ok=false;
        
        if (ok) 
        {
            double lgt,lat,alt;
            ECEF_to_lgt_lat_alt(ecef[w*3*y+3*x], 
                                ecef[w*3*y+3*x+1], 
                                ecef[w*3*y+3*x+2],
                                &lgt,&lat,&alt);
            get_utm_coord(utm_coord, lat, lgt, alt, zone);
            fwrite( utm_coord, sizeof(double), point_size, ply_file);
        }
    }

    free(utm_coord);
    free(ecef);
    fclose(ply_file);
    return 0;
}
