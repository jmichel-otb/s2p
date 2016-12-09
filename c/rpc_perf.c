#include "stdio.h"
#include "stdlib.h"
#include "rpc.h"
#include "refine_rpc.h"

int main_rpc_refiner(int c, char *v[])
{
    
    if (c != 4) {
        fprintf(stderr, "usage:\n\t"
                "%s rpc.xml rpc_refined.xml tie_points.txt "
              // 0     1         2              3
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }
    
    // parse arguments
    struct rpc rpc_coef;
    read_rpc_file_xml(&rpc_coef,v[1]);
    
    struct rpc refined_rpc_coef;
    read_rpc_file_xml(&refined_rpc_coef,v[2]);
    
    unsigned int nb_tie_points;
    get_nb_tie_points(v[3], &nb_tie_points);
    Tie_point* list_tie_points = (Tie_point*) malloc(nb_tie_points*sizeof(Tie_point));
    get_tie_points(v[3], list_tie_points, nb_tie_points);
       
   // Perfs
    printf("Perf direct model\n");
    double perf_before = perf_rpc(&rpc_coef,list_tie_points, nb_tie_points);
    double perf_after = perf_rpc(&refined_rpc_coef,list_tie_points, nb_tie_points);
    printf("perf (before)= %f\n",perf_before);
    printf("perf (after)= %f\n",perf_after);

    printf("Perf indirect model\n");
    double perfi_before = perf_rpci(&rpc_coef,list_tie_points, nb_tie_points);
    double perfi_after = perf_rpci(&refined_rpc_coef,list_tie_points, nb_tie_points);
    printf("perfi (before)= %f\n",perfi_before);
    printf("perfi (after)= %f\n",perfi_after);

    return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
    return main_rpc_refiner(c, v);
}
