#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"


int main_nan_generators(int c, char *v[])
{
    if (c != 4) {
        fprintf(stderr, "usage:\n\t"
                "%s out_img_full_path, tw, th"
              // 0         1         2  3   
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }
    
    char *out_img_full_path=v[1];
    int width=atoi(v[2]);
    int height=atoi(v[3]);
    
    // build outputs paths and alloc mem
    //char fout[1000];
    //sprintf(fout,"%s/nan.tif",out_img_full_path);
    float *nanimg = (float *) malloc(3*width*height*sizeof(float));
    
    // nans
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            for(int t=0; t<3; t++)
                    nanimg[width*3*y+3*x+t] = NAN;
        
    // save 
    iio_save_image_float_vec(out_img_full_path, nanimg, width, height,3);
    
    // clean mem
    free(nanimg);
    
    return 0;
}

int main(int c, char *v[])
{
    return main_nan_generators(c, v);
}
