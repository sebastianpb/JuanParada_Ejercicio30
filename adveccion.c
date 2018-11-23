#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>


void initial(double delta_x, long double *x, long double *u, double L) {
    double dN=(L/delta_x) + 1;
    int N=(int)dN;
    double fr=L/(N-1);
    int i;

    for(i = 0; i < N; i++)
    {
        x[i]=(double)i*fr;
        if(x[i]<2)
        {
            u[i]=1;
        }
        else{
            u[i]=0;
        }
    }
}


void flux(long double *u, long double *F, int N) {
    int i;

    for(i = 0; i < N; i++)
    {
        F[i] = 0.5*pow(u[i],2);
    }
}

    
void Lax(long double *u, long double *u_final, double t_max, double delta_t, double delta_x, int N) {
    double fr=t_max/delta_t;
    int N_t=(int)fr;
    int i,j;
    long double F[N];    

    for(i = 0; i < N_t; i++)
    {
        flux(u,F,N);
        
        for(j=0;j<N-2;j++){
            u_final[j+1]=0.5*(u[j+2]+u[j]);
            u_final[j+1]-=(0.5*delta_t/delta_x)*(F[j+2] - F[j]) ;
        }
        
        for(j = 0; j < N; j++){
            u[j]=u_final[j];
        }
    }
}

int main(int argc, char **argv){

    int i,j;
    double L=4.0;
    double delta_x=0.05;
    double dN=(L/delta_x) + 1;
    int N=(int)dN;
    
    long double u[N];
    long double x[N];
    long double u_final[N];
    
    double t_max_values[4]={0, 0.5, 1, 2.0};
    double delta_t=0.5*delta_x;
    
    omp_set_num_threads(4);
    #pragma omp parallel
    {    
        
        int thread_id=omp_get_thread_num();
        int thread_count=omp_get_num_threads();
        printf("Hello from thread number: %d out of: %d\n",thread_id,thread_count);
        
        FILE *out;
        double t_max=t_max_values[thread_id];
        
            // Inicializar todo en cero
        for(i=0;i<N;i++){
            u[i]=0;
            x[i]=0;
            u_final[N]=0;
        }    

        initial(delta_x, x,u,L); 

        for(i=0;i<N;i++){
            u_final[i]=u[i];
        }

        Lax(u, u_final,  t_max, delta_t, delta_x, N);


        // Imprimir resultados

        /*
        for(i=0;i<N;i++)
        {
            printf("u[%d]: %Lf\n",i,u[i]);
        }
        for(i=0;i<N;i++)
        {
            printf("x[%d]= %.2Lf\n",i,x[i]);
        }
        */
        for(i=0;i<N;i++)
        {
            printf("u_final[%d]= %.8Lf\n",i,u_final[i]);
        }

        // Save into dat file
        char filename[128];
        double num=t_max*10;
        sprintf(filename, "adveccion_%d.dat", (int)num);

        if(!(out = fopen(filename, "w"))){
            fprintf(stderr, "Problema abriendo el archivo\n");
            exit(1);
        }
        for(i=0;i<N;i++){
            fprintf(out, "%Lf %Lf\n", u_final[i],x[i]);
        }
        fclose(out);      
        
    }
    
     // Save x vals into dat file
    char filename[128];
    FILE *out;
    sprintf(filename, "adveccion_x.dat");

    if(!(out = fopen(filename, "w"))){
        fprintf(stderr, "Problema abriendo el archivo\n");
        exit(1);
    }
    for(i=0;i<N;i++){
        fprintf(out, "%Lf\n", x[i]);
    }
    fclose(out); 
    
}