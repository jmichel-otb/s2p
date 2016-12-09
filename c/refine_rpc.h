// rational polynomial coefficient stuff
#ifndef _REFINE_RPC_H
#define _REFINE_RPC_H

typedef struct {
    double x;
    double y;
    double lgt;
    double lat;
    double alt;
} Tie_point;

int get_nb_tie_points(char *filename, unsigned int *nb_tie_points);

int get_tie_points(char *filename, Tie_point* list_tie_points, unsigned int nb_tie_points);

double perf_rpci(struct rpc *rpc_coef,Tie_point* list_tie_points, unsigned int nb_tie_points);

double perf_rpc(struct rpc *rpc_coef,Tie_point* list_tie_points, unsigned int nb_tie_points);

double * get_address(struct rpc *rpc_coef,int i);

double * get_addressi(struct rpc *rpc_coef,int i);

double perf_deriv(double **addr, int i, double step, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points);

void gradient(double **addr, int size, double h, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points,
double *out);

void update(double **addr, int size, double *gradient,double step);

void gradient_descent(double **addr, int size, double step_deriv, double step_grad, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points, int nb_iter);

int write_rpc_coef(char *filename,double **addr);

#endif // _REFINE_RPC_H
