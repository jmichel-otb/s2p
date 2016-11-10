#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "rpc.h"
#include "refine_rpc.h"


int get_tie_points(char *filename, Tie_point* list_tie_points, unsigned int nb_tie_points)
{
    double a,b,c,d,e;
    FILE *fic=fopen(filename,"r");
    if (fic)
    {
        if (nb_tie_points>0)
        {
            for(unsigned int t=0; t<nb_tie_points; t++)
                fscanf(fic,"%lf %lf %lf %lf %lf\n",&list_tie_points[t].x,
                            &list_tie_points[t].y,
                            &list_tie_points[t].lgt,
                            &list_tie_points[t].lat,
                            &list_tie_points[t].alt);
            fclose(fic);
            return 0;
        }
    }
    else
        return 1;
}

int get_nb_tie_points(char *filename, unsigned int *nb_tie_points)
{
    *nb_tie_points = 0;
    FILE *fic=fopen(filename,"r");
    char chartab[1000];
    if (fic)
    {
        while(fgets(chartab,1000,fic))
            (*nb_tie_points)++;
        fclose(fic); 
        return 0;
    }
    else
        return 1;
}

double perf_rpci(struct rpc *rpc_coef,Tie_point* list_tie_points, unsigned int nb_tie_points)
{
    double pos[2],diffx,diffy;
    double diffx_moy=0.,diffy_moy=0.,tot_moy=0.0;
    for(unsigned int t=0;t<nb_tie_points;t++)
    {
        eval_rpci(pos, rpc_coef, list_tie_points[t].lgt, list_tie_points[t].lat, list_tie_points[t].alt);
        diffx = pow(pos[0]-list_tie_points[t].x,2.0);
        diffy = pow(pos[1]-list_tie_points[t].y,2.0);
        //diffx_moy += diffx;
        //diffy_moy += diffy;
        tot_moy += diffx + diffy;
        //printf("%f %f --> %f %f  (%f %f)\n",list_tie_points[t].x,list_tie_points[t].y,pos[0],pos[1],diffx,diffy);
    }
    //diffx_moy = sqrt( diffx_moy / ( (double) nb_tie_points) );
    //diffy_moy = sqrt( diffy_moy / ( (double) nb_tie_points) );
    tot_moy = sqrt( tot_moy / ( (double) nb_tie_points) );
    //printf("RMS.x = %f  RMS.y = %f RMS.tot = %f \n",diffx_moy,diffy_moy,tot_moy);
    return tot_moy;
}


double perf_rpc(struct rpc *rpc_coef,Tie_point* list_tie_points, unsigned int nb_tie_points)
{
    double pos[2],difflgt,difflat;
    double difflgt_moy=0.,difflat_moy=0.,tot_moy=0.0;
    for(unsigned int t=0;t<nb_tie_points;t++)
    {
        eval_rpc(pos, rpc_coef, list_tie_points[t].lgt, list_tie_points[t].lat, list_tie_points[t].alt);
        difflgt = pow(pos[0]-list_tie_points[t].lgt,2.0);
        difflat = pow(pos[1]-list_tie_points[t].lat,2.0);
        //difflgt_moy += difflgt;
        //difflat_moy += difflat;
        tot_moy += difflgt + difflat;
        //printf("%f %f --> %f %f  (%f %f)\n",list_tie_points[t].lgt,list_tie_points[t].lat,pos[0],pos[1],difflgt,difflat);
    }
    //difflgt_moy = sqrt( difflgt_moy / ( (double) nb_tie_points) );
    //difflat_moy = sqrt( difflat_moy / ( (double) nb_tie_points) );
    tot_moy = sqrt( tot_moy / ( (double) nb_tie_points) );
    //printf("RMS.lgt = %f  RMS.lat = %f RMS.tot = %f \n",difflgt_moy,difflat_moy,tot_moy);
    return tot_moy;
}


double * get_address(struct rpc *p, int i)
{
    if (false);
    // pleiades tags
	else if (i==0)   return &p->numx[0];
	else if (i==1)   return &p->numx[1];
	else if (i==2)   return &p->numx[2];
	else if (i==3)   return &p->numx[3];
	else if (i==4)   return &p->numx[4];
	else if (i==5)   return &p->numx[5];
	else if (i==6)   return &p->numx[6];
	else if (i==7)   return &p->numx[7];
	else if (i==8)   return &p->numx[8];
	else if (i==9)   return &p->numx[9];
	else if (i==10)   return &p->numx[10];
	else if (i==11)   return &p->numx[11];
	else if (i==12)   return &p->numx[12];
	else if (i==13)   return &p->numx[13];
	else if (i==14)   return &p->numx[14];
	else if (i==15)   return &p->numx[15];
	else if (i==16)   return &p->numx[16];
	else if (i==17)   return &p->numx[17];
	else if (i==18)   return &p->numx[18];
	else if (i==19)   return &p->numx[19];
    
    else if (i==20)   return &p->denx[0];
	else if (i==21)   return &p->denx[1];
	else if (i==22)   return &p->denx[2];
	else if (i==23)   return &p->denx[3];
	else if (i==24)   return &p->denx[4];
	else if (i==25)   return &p->denx[5];
	else if (i==26)   return &p->denx[6];
	else if (i==27)   return &p->denx[7];
	else if (i==28)   return &p->denx[8];
	else if (i==29)   return &p->denx[9];
	else if (i==30)   return &p->denx[10];
	else if (i==31)   return &p->denx[11];
	else if (i==32)   return &p->denx[12];
	else if (i==33)   return &p->denx[13];
	else if (i==34)   return &p->denx[14];
	else if (i==35)   return &p->denx[15];
	else if (i==36)   return &p->denx[16];
	else if (i==37)   return &p->denx[17];
	else if (i==38)   return &p->denx[18];
	else if (i==39)   return &p->denx[19];
    
    else if (i==40)   return &p->numy[0];
	else if (i==41)   return &p->numy[1];
	else if (i==42)   return &p->numy[2];
	else if (i==43)   return &p->numy[3];
	else if (i==44)   return &p->numy[4];
	else if (i==45)   return &p->numy[5];
	else if (i==46)   return &p->numy[6];
	else if (i==47)   return &p->numy[7];
	else if (i==48)   return &p->numy[8];
	else if (i==49)   return &p->numy[9];
	else if (i==50)   return &p->numy[10];
	else if (i==51)   return &p->numy[11];
	else if (i==52)   return &p->numy[12];
	else if (i==53)   return &p->numy[13];
	else if (i==54)   return &p->numy[14];
	else if (i==55)   return &p->numy[15];
	else if (i==56)   return &p->numy[16];
	else if (i==57)   return &p->numy[17];
	else if (i==58)   return &p->numy[18];
	else if (i==59)   return &p->numy[19];
    
    else if (i==60)   return &p->deny[0];
	else if (i==61)   return &p->deny[1];
	else if (i==62)   return &p->deny[2];
	else if (i==63)   return &p->deny[3];
	else if (i==64)   return &p->deny[4];
	else if (i==65)   return &p->deny[5];
	else if (i==66)   return &p->deny[6];
	else if (i==67)   return &p->deny[7];
	else if (i==68)   return &p->deny[8];
	else if (i==69)   return &p->deny[9];
	else if (i==70)   return &p->deny[10];
	else if (i==71)   return &p->deny[11];
	else if (i==72)   return &p->deny[12];
	else if (i==73)   return &p->deny[13];
	else if (i==74)   return &p->deny[14];
	else if (i==75)   return &p->deny[15];
	else if (i==76)   return &p->deny[16];
	else if (i==77)   return &p->deny[17];
	else if (i==78)   return &p->deny[18];
	else if (i==79)   return &p->deny[19];
    
    else if (i==80)   return &p->scale[0];
	else if (i==81)   return &p->scale[1];
	else if (i==82)   return &p->scale[2];
	else if (i==83)   return &p->offset[0];
	else if (i==84)   return &p->offset[1];
	else if (i==85)   return &p->offset[2];
    
    else if (i==86)   return &p->iscale[0];
	else if (i==87)   return &p->iscale[1];
	else if (i==88)   return &p->ioffset[0];
	else if (i==89)   return &p->ioffset[1];
}

double * get_addressi(struct rpc *p, int i)
{
    if (false);
    // pleiades tags
	else if (i==0)   return &p->inumx[0];
	else if (i==1)   return &p->inumx[1];
	else if (i==2)   return &p->inumx[2];
	else if (i==3)   return &p->inumx[3];
	else if (i==4)   return &p->inumx[4];
	else if (i==5)   return &p->inumx[5];
	else if (i==6)   return &p->inumx[6];
	else if (i==7)   return &p->inumx[7];
	else if (i==8)   return &p->inumx[8];
	else if (i==9)   return &p->inumx[9];
	else if (i==10)   return &p->inumx[10];
	else if (i==11)   return &p->inumx[11];
	else if (i==12)   return &p->inumx[12];
	else if (i==13)   return &p->inumx[13];
	else if (i==14)   return &p->inumx[14];
	else if (i==15)   return &p->inumx[15];
	else if (i==16)   return &p->inumx[16];
	else if (i==17)   return &p->inumx[17];
	else if (i==18)   return &p->inumx[18];
	else if (i==19)   return &p->inumx[19];
    
    else if (i==20)   return &p->idenx[0];
	else if (i==21)   return &p->idenx[1];
	else if (i==22)   return &p->idenx[2];
	else if (i==23)   return &p->idenx[3];
	else if (i==24)   return &p->idenx[4];
	else if (i==25)   return &p->idenx[5];
	else if (i==26)   return &p->idenx[6];
	else if (i==27)   return &p->idenx[7];
	else if (i==28)   return &p->idenx[8];
	else if (i==29)   return &p->idenx[9];
	else if (i==30)   return &p->idenx[10];
	else if (i==31)   return &p->idenx[11];
	else if (i==32)   return &p->idenx[12];
	else if (i==33)   return &p->idenx[13];
	else if (i==34)   return &p->idenx[14];
	else if (i==35)   return &p->idenx[15];
	else if (i==36)   return &p->idenx[16];
	else if (i==37)   return &p->idenx[17];
	else if (i==38)   return &p->idenx[18];
	else if (i==39)   return &p->idenx[19];
    
    else if (i==40)   return &p->inumy[0];
	else if (i==41)   return &p->inumy[1];
	else if (i==42)   return &p->inumy[2];
	else if (i==43)   return &p->inumy[3];
	else if (i==44)   return &p->inumy[4];
	else if (i==45)   return &p->inumy[5];
	else if (i==46)   return &p->inumy[6];
	else if (i==47)   return &p->inumy[7];
	else if (i==48)   return &p->inumy[8];
	else if (i==49)   return &p->inumy[9];
	else if (i==50)   return &p->inumy[10];
	else if (i==51)   return &p->inumy[11];
	else if (i==52)   return &p->inumy[12];
	else if (i==53)   return &p->inumy[13];
	else if (i==54)   return &p->inumy[14];
	else if (i==55)   return &p->inumy[15];
	else if (i==56)   return &p->inumy[16];
	else if (i==57)   return &p->inumy[17];
	else if (i==58)   return &p->inumy[18];
	else if (i==59)   return &p->inumy[19];
    
    else if (i==60)   return &p->ideny[0];
	else if (i==61)   return &p->ideny[1];
	else if (i==62)   return &p->ideny[2];
	else if (i==63)   return &p->ideny[3];
	else if (i==64)   return &p->ideny[4];
	else if (i==65)   return &p->ideny[5];
	else if (i==66)   return &p->ideny[6];
	else if (i==67)   return &p->ideny[7];
	else if (i==68)   return &p->ideny[8];
	else if (i==69)   return &p->ideny[9];
	else if (i==70)   return &p->ideny[10];
	else if (i==71)   return &p->ideny[11];
	else if (i==72)   return &p->ideny[12];
	else if (i==73)   return &p->ideny[13];
	else if (i==74)   return &p->ideny[14];
	else if (i==75)   return &p->ideny[15];
	else if (i==76)   return &p->ideny[16];
	else if (i==77)   return &p->ideny[17];
	else if (i==78)   return &p->ideny[18];
	else if (i==79)   return &p->ideny[19];
    
    else if (i==80)   return &p->iscale[0];
	else if (i==81)   return &p->iscale[1];
	else if (i==82)   return &p->iscale[2];
	else if (i==83)   return &p->ioffset[0];
	else if (i==84)   return &p->ioffset[1];
	else if (i==85)   return &p->ioffset[2];
    
    else if (i==86)   return &p->scale[0];
	else if (i==87)   return &p->scale[1];
	else if (i==88)   return &p->offset[0];
	else if (i==89)   return &p->offset[1];
}

double perf_deriv(double **addr, int i, double h, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points)
{
    double val1,val2;
    double old = *addr[i];
    
    // f(...,xi+h,...)
    *addr[i] = old + h;
    if (direct)
        val1 = perf_rpc(rpc_coef,list_tie_points, nb_tie_points);
    else
        val1 = perf_rpci(rpc_coef,list_tie_points, nb_tie_points);
    
    // f(...,x-h,...)    
    *addr[i] = old - h;
    if (direct)
        val2 = perf_rpc(rpc_coef,list_tie_points, nb_tie_points);
    else
        val2 = perf_rpci(rpc_coef,list_tie_points, nb_tie_points);
    
    // back to original value    
    *addr[i] = old;
    
    // [ f(...,xi+h,...) - f(...,x-h,...) ] / 2h
    return (val1-val2)/(2.0*h);
}

void gradient(double **addr, double h, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points,
double *out)
{
    for(int i=0;i<90;i++)
        out[i] = perf_deriv(addr, i, h, direct,
        rpc_coef, list_tie_points, nb_tie_points);
    
}

double norm_gradient(double *gradient)
{
    double sum=0;
    for(int i=0;i<90;i++)
        sum += pow(gradient[i],2.0);
    return sqrt(sum);
}

void update(double **addr, double *gradient,double step)
{
    for(int i=0;i<90;i++)
        *addr[i] += step*gradient[i];
}


void gradient_descent(double **addr, double step_deriv, double step_grad, bool direct,
struct rpc *rpc_coef, Tie_point* list_tie_points, unsigned int nb_tie_points, int nb_iter)
{
    double gradient_val[90];
    double norm,perf;
    
    for(int i=0;i<nb_iter;i++)
    {
        gradient(addr,step_deriv, direct,
            rpc_coef, list_tie_points, nb_tie_points, gradient_val);
        norm = norm_gradient(gradient_val);
        update(addr, gradient_val,step_grad);
        
        if (direct)
            perf = perf_rpc(rpc_coef,list_tie_points, nb_tie_points);
        else
            perf = perf_rpci(rpc_coef,list_tie_points, nb_tie_points);
        printf("perf = %f  norm gradient = %f  (i=%d)\n",perf,norm,i);
    }
}

int write_rpc_coef(char *filename,double **addr)
{
    FILE *fic = fopen(filename,"w");
    if (fic)
    {
        for(int i=0;i<90;i++)
            fprintf(fic,"%.19f\n",*addr[i]);
        fclose(fic);
        return 0;
    }
    else
        return 1;
}
