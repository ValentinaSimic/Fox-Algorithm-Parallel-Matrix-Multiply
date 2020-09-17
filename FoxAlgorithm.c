#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>



typedef struct {
    MPI_Comm grid_comm; 
    MPI_Comm row_comm;  
    MPI_Comm col_comm;  
    int p;         
    int grid_dim;       
    int my_row;        
    int my_col;         
    int my_rank;       
} Grid;



void grid_init(Grid *grid)
{
    int old_rank;
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int remain_dims[2];

   
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

    grid->grid_dim = (int)pow(grid->p,0.5);
    
    if (grid->grid_dim * grid->grid_dim != grid->p) {
        printf("Dimenzija resetke mora biti kvadratni koren od ukupnog broja procesa!\n");
        exit(-1);
    }
  
    dimensions[0] = grid->grid_dim;
    dimensions[1] = grid->grid_dim;
    wrap_around[0] =1;
    wrap_around[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid->grid_comm));
   
    MPI_Comm_rank(grid->grid_comm, &(grid->my_rank));
  
    MPI_Cart_coords(grid->grid_comm, grid->my_rank, 2, coordinates);
 
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

   
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(grid->grid_comm, remain_dims, &(grid->row_comm));

  
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(grid->grid_comm, remain_dims, &(grid->col_comm));
}


void matrix_multiply(double *A, double *B, double *C, int n)
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int z=0; z<n; z++) {
                C[i*n+j] += A[i*n+z] * B[z*n+j];
            }
        }
    }
}


void matrix_print(double *A, int n)
{
    printf("***************************\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%5d ", A[i*n+j]);
        }
        printf("\n");
    }
      printf("***************************\n\n");
}




void FoxAlgorithm(double *A, double *B, double *C, int n, Grid *grid)
{
    double *buff_A = (double*)calloc(n * n, sizeof(double));
    MPI_Status status;
    int root;
    int source = (grid->my_row + 1) % grid->grid_dim;
    int destination = (grid->my_row -1 + grid->grid_dim) % grid->grid_dim;

    for (int stage = 0; stage < grid->grid_dim; ++stage) {
        root = (grid->my_row + stage) % grid->grid_dim;
        if (root == grid->my_col) {
            MPI_Bcast(A, n*n, MPI_DOUBLE, root, grid->row_comm);
            matrix_multiply(A, B, C, n);
        } else {
            MPI_Bcast(buff_A, n*n, MPI_DOUBLE, root, grid->row_comm);
            matrix_multiply(buff_A, B, C, n);
        }
        MPI_Sendrecv_replace(B, n*n, MPI_DOUBLE, destination, 0, source, 0, grid->col_comm, &status);
    }
}



int main(int argc, char **argv) {

    if (argc <1) {
        printf("Morate zadati dimenziju matrica kao argument komandne linije!\n");
        return -1;
    }
    double *A, *B, *result;

    int n=(int)atoi(argv[1]);
    
    MPI_Init(&argc, &argv);
   
    Grid grid;
    grid_init(&grid);
   
    if (n % grid.grid_dim != 0) {
        printf("Kod Fox algoritma dimenzija matrice mora biti deljiva sa kvadratnim korenom broja procesa !\n");
        exit(-1);
    }
   
    if (grid.my_rank == 0) {
    
        A = (double *)malloc(n*n * sizeof(double));
        B = (double *)malloc(n*n * sizeof(double));
        result = (double *)calloc(n*n, sizeof(double));
 
       //printf("Unesite prvu matricu dimenzije %dx%d \n", n, n);
       for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            //scanf("%d",&A[i*n+j]);
            A[i*n+j]= (double) (rand() % (5 - 1 + 1)) + 1; 
        }
       }
       // printf("Unesite drugu matricu dimenzije %dx%d\n", n,n);
       for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            //scanf("%d",&B[i*n+j]);
            B[i*n+j]=(double) (rand() % (5 - 1 + 1)) + 1;
        }
       }
       // matrix_print(A,matrix_size);
       //matrix_print(B, matrix_size);
    }
    int local_n = n / grid.grid_dim;
    double local_A[local_n*local_n];
    double local_B[local_n*local_n];
    double *local_result = (double*)calloc(local_n*local_n, sizeof(double));
    double A_temp_B[n*n];
   
   
    if (grid.my_rank == 0) { 

        int i = 0;
        int j = 0;
        for(int row_num = 0; row_num < grid.grid_dim; row_num++){
          for(int col_num= 0; col_num < grid.grid_dim; col_num++){
            j = col_num*local_n + row_num*n*local_n;
            for(int k = 0; k < local_n*local_n; k++){
              if(k% local_n== 0 && k != 0){
                j = j + (n-local_n);
              }

              A_temp_B[i] = A[j];
            
              i++;
              j++;   
            }
          }
        }

        i = 0;
        j = 0;
        for(int row_num = 0; row_num < grid.grid_dim; row_num++){
          for(int col_num= 0; col_num < grid.grid_dim; col_num++){
            j = col_num*local_n + row_num*n*local_n;
            for(int k = 0; k < local_n*local_n; k++){
              if(k% local_n == 0 && k != 0){
                j = j + (n-local_n);
              }

              A[i] = B[j];
            
              i++;
              j++;
            }
          }
        }
     //matrix_print(A_temp_B, matrix_size);
     //matrix_print(A, matrix_size);
    
     }
   MPI_Scatter(A_temp_B, (local_n*local_n), MPI_DOUBLE, &local_A, (local_n*local_n), MPI_DOUBLE,0, MPI_COMM_WORLD);
   MPI_Scatter(A, local_n*local_n, MPI_DOUBLE, &local_B, local_n*local_n, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    

 
    double start_time, finish_time;

    MPI_Barrier(grid.grid_comm);
    if (grid.my_rank == 0) {
        start_time = MPI_Wtime();
    }
 
    FoxAlgorithm(local_A, local_B, local_result, local_n, &grid);
    
    MPI_Barrier(grid.grid_comm);
    if (grid.my_rank == 0) {
        finish_time = MPI_Wtime() - start_time;
    }
   
   MPI_Gather(local_result, local_n*local_n, MPI_DOUBLE, result, local_n*local_n,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    if (grid.my_rank == 0) {
       
    
    //matrix_print(result, matrix_size);
    printf("Vreme izvrsavanja: %.10lf, broj procesa: %d, dimenzija matrica: %d\n", finish_time, grid.p, n);
    free(A);
    free(B);
    free(result);
        
    }
    
    free(local_result);
    MPI_Finalize();
    
    return 0;
}