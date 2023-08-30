#include<iostream>
#include<mpi.h>
#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/core.h>

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    //obtener el rank y numero de procesos
    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //rank
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //procesos

    int M;
    float *b = new float[M];
    float *A = new float[M * M];
    float *c = new float[M];


    if (rank == 0) {
        M = 32;
//Inicializar
        for (int i = 0; i < M * M; i++) {
            A[i] = i;
        }
        for (int i = 0; i < M; i++) {
            b[i] = i;
            c[i] = 0;
        }
//Multiplicacion Serial-----------------------------------------------------------------

        /* M=3
         * A[0,1,2,3,4,5,6,7,8]
         * b[0,1,2]
         * c[(0+1+4),(0+4+10),(0+7+16)] = c[5,14,23]
         */

        for (int j = 0; j < M * M; j += M) {
            float tmp = 0;
            for (int i = 0; i < M; i++) {
                tmp += A[i + j] * b[i];
            }
            c[j / M] = tmp;
        }
//---------------------------------------------------------------------------------------


    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b[0], M, MPI_FLOAT, 0, MPI_COMM_WORLD);


    if (rank == 0) {

        /* M=32
         * nprocs=4
         * nfilas=8
         * A[32*32]
         * rank0= A[0-7 y 8-15]=0-15
         * ran1= A[16-31] (nrank)*M*nfilas =16
         * rank2=A[32-47] (2*8*2)=32
         * rank3=A[48-63] (3*8*2)=48
         * bloque= M*nfilas
         */
        int nfilas = M / nprocs;
        int bloque = M * nfilas;
        for (int nRank = 1; nRank < nprocs; nRank++) {
            MPI_Send(&A[nRank * M * nfilas], bloque, MPI_FLOAT, nRank, 0, MPI_COMM_WORLD);
        }

        for (int j = 0; j < M * nfilas; j += M) {
            float tmp = 0;
            for (int i = 0; i < M; i++) {
                tmp += A[i + j] * b[i];
            }
            c[j / M] = tmp;
        }

        for (int nRank = 1; nRank < nprocs; nRank++) {
            MPI_Recv(&c[nRank * nfilas], nfilas, MPI_FLOAT, nRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }




    } else {

        int nfilas = M / nprocs;
        int bloque = M * nfilas;

        float *ATemp = new float[bloque];
        MPI_Recv(&ATemp[0], bloque, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        float *cTmp = new float[nfilas];

        for (int j = 0; j < M * nfilas; j += M) {
            float tmp = 0;
            for (int i = 0; i < M; i++) {
                tmp += ATemp[i + j] * b[i];
            }
            cTmp[j / M] = tmp;
        }

        MPI_Send(&cTmp[0], nfilas, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);


    }

    MPI_Finalize();
    return 0;
}