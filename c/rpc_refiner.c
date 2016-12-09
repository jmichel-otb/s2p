#include "stdio.h"
#include "stdlib.h"
#include "rpc.h"
#include "refine_rpc.h"

int main_rpc_refiner(int c, char *v[])
{
    
    if (c != 9) {
        fprintf(stderr, "usage:\n\t"
                "%s rpc.xml tie_points.txt step_deriv step_grad nb_iter bool_direct bool_first_80_coefs out.txt"
              // 0     1      2               3           4       5          6               7             8
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }
    
    // parse arguments
    struct rpc rpc_coef;
    read_rpc_file_xml(&rpc_coef,v[1]);
    
    unsigned int nb_tie_points;
    get_nb_tie_points(v[2], &nb_tie_points);
    Tie_point* list_tie_points = (Tie_point*) malloc(nb_tie_points*sizeof(Tie_point));
    get_tie_points(v[2], list_tie_points, nb_tie_points);
    
    double step_deriv = atof(v[3]);
    double step_grad = atof(v[4]);
    int nb_iter = atoi(v[5]);
    bool direct = (bool) atoi(v[6]);
    bool only_first_80_coefs = (bool) atoi(v[7]);
    char * fname_out = v[8];
    
    // convert struct rpc to a tab of 90 double
    double *address[90];
    double *addressi[90];
    for(int i=0;i<90;i++)
    {
        address[i] = get_address(&rpc_coef,i);
        addressi[i] = get_addressi(&rpc_coef,i);
    }
    
    int size=90;
    if (only_first_80_coefs)
        size=80;
    
    // refining
    if ( direct)
    {   
        printf("Refining direct model\n");
        double perf_before = perf_rpc(&rpc_coef,list_tie_points, nb_tie_points);
        gradient_descent(address, size, step_deriv, step_grad, true,
            &rpc_coef, list_tie_points, nb_tie_points, nb_iter);
        double perf_after = perf_rpc(&rpc_coef,list_tie_points, nb_tie_points);
        printf("perf (before)= %f\n",perf_before);
        printf("perf (after)= %f\n",perf_after);
        write_rpc_coef(fname_out,address);
    }
    else
    {
        printf("Refining indirect model\n");
        double perfi_before = perf_rpci(&rpc_coef,list_tie_points, nb_tie_points);
        gradient_descent(addressi, size, step_deriv, step_grad, false,
            &rpc_coef, list_tie_points, nb_tie_points, nb_iter);
        double perfi_after = perf_rpci(&rpc_coef,list_tie_points, nb_tie_points);
        printf("perfi (before)= %f\n",perfi_before);
        printf("perfi (after)= %f\n",perfi_after);
        write_rpc_coef(fname_out,addressi);
    }

    return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
    return main_rpc_refiner(c, v);
}
