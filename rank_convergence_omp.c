/* 	program dot3mpi.c */
/* Illustrates dot product via mpi_gather.*/
#include <stdio.h>
#include <omp.h>
/* Includes the mpi C library. */
#include "mpi.h"
#include <limits.h>
#include <stdlib.h>

/* contants defined
* n : size of sequence of observations
* m : number of hidden states
* t : size of the current observation
* A : Transmission matrix assumed to be time invariant (m*m)
* V : The Viterbi output matrix of size (m*t)
* B : Emission probablilty matrix (m*n)
* init : Initial probability matrix(m)
*/

struct interm_type{ 
        double val; 
        int   index; 
};

      int main(int argc, char* argv[]) {
      int n=100,m=100,t=100;
      
      
      int  my_rank,p,source,dest,loc_n;
      int  i,en,bn;
      MPI_Status  status;
      MPI_Comm comm = MPI_COMM_WORLD;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(comm,&my_rank);
      MPI_Comm_size(comm,&p);

      int root=0;
      
/* Initializes mpi, gets the rank of the processor, my_rank, */
/* and number of processors, p. */
      
      if(my_rank==0)
      scanf("%d %d %d",&n,&m,&t);
      MPI_Bcast(&n,1,MPI_INT,root,comm);
      MPI_Bcast(&m,1,MPI_INT,root,comm);
      MPI_Bcast(&t,1,MPI_INT,root,comm);

      double A[m][m], B[m][n], init[m];
      int result[t];
      int argmnt[m][t];
      int observe[t];

/* Each processor computes a local dot product */
      loc_n = m/p;
      int send_count[p]; // Array used to store the count of data send to each processor
      int display[p]; // array specifies the displacement from which to take the outgoing data to process
      display[0] = 0;
      for (int i = 0; i < p; ++i)
      {
         send_count[i] = loc_n;
         if (i < m%p)
         {
            send_count[i]++;
         }
         if(i<p-1)
            display[i+1] = display[i] + send_count[i];
      }

      bn = (my_rank)*loc_n;
      en = bn + loc_n-1;
      double proc_array[loc_n+2]; // Local array for each processo to receive the scattered data
      struct interm_type interm_prod[m]; // array to store the maximum value of the loc_n states along with the correesponding index
      struct interm_type red[m]; // array used to receive the global maximum along with indices for each state, present in the root process
      double sol[m]; // vector used to store the current stage
      printf("my_rank = %d loc_n = %d\n",my_rank,loc_n);
      printf("my_rank = %d bn = %d\n",my_rank,bn);
      printf("my_rank = %d en = %d\n",my_rank,en);
      
      if(my_rank==0)
      {
          /* Taking input */
          for(int i=0;i<m;i++)
          {
              for(int j=0;j<m;j++)
              {
                scanf("%lf",&A[i][j]);
              }
          }
          
          for(int i=0;i<m;i++){
            scanf("%lf",&init[i]);
          }
        //init[i]=((double)(rand()%20))/100;
          for(int i=0;i<m;i++)
          {
              for(int j=0;j<n;j++)
              {
                scanf("%lf",&B[i][j]);
              }
          } 
        //Input the observation vector(must have numbers between 0...n-1)   
          for(int i=0;i<t;i++)
            scanf("%d",&observe[i]);
          for(int i=0;i<m;i++)
            sol[i] = init[i]*B[i][observe[0]];
            //MPI_Barrier(comm); // Initialize the first stage 
          for(int i=0;i<p;i++)
          {
            printf("%d %d\n",send_count[i],display[i]);
          }


    }
    MPI_Bcast(A,m*m,MPI_DOUBLE,root,comm);
    

    // /*if(my_rank==1)
    //   {
    //       /* Taking input */
    //       for(int i=0;i<m;i++)
    //       {
    //           for(int j=0;j<m;j++)
    //           {
    //             printf("%lf ",A[i][j]);
    //           }
    //           printf("\n");
    //       }
    //   }

      for(int i=1;i<t;i++)
      {
          //Scattering send_count[i] elements of the sol vector to the processors 
          // if(my_rank==0)
          // {
          //   for(int k=0;k<m;k++)
          //     printf("%lf ",sol[k]);
          //   printf("\n");
          // }
        

          MPI_Scatterv(sol,send_count,display,MPI_DOUBLE,proc_array,loc_n+2,MPI_DOUBLE,0,comm);
          //#pragma omp parallel for
          for(int k=0;k<m;k++)
          {
      			  double mx = INT_MIN;
      			  int indx=0;
      			  for (int j = 0; j < send_count[my_rank]; j++) // // collect local_max
      			  {
                if(i==1 && my_rank==0){
                  printf("%d %d %d %lf %lf \n",j,display[my_rank]+j,k,proc_array[j],proc_array[j]*A[display[my_rank]+j][k]);
                }
      				  if(proc_array[j]*A[display[my_rank]+j][k] > mx)
                {
          				  mx = proc_array[j]*A[display[my_rank]+j][k];
          				  indx = display[my_rank] + j;
      				  }
      			  }
              if(my_rank==0 && i==1)
                printf("\n");

      			  interm_prod[k].val = mx;
      			  interm_prod[k].index = indx;
		      }


		  //After finding the local max and the corresponding index perform reduction 
		  MPI_Reduce( interm_prod, red, m, MPI_DOUBLE_INT, MPI_MAXLOC, root, comm );
		  
		  //Copy the values to the pred array and to the sol vector which will be scattered
		  if(my_rank==0)
		  {
  			for(int k=0;k<m;k++)
  			{
    				sol[k] = red[k].val * B[k][i];
    				argmnt[k][i] = red[k].index;
  			}
  		}
    }
    // Now the argument matrix has been populated 
	

  if(my_rank ==0)
	{
    	double mx = INT_MIN;
      int indx = 0;
      for(int i=0;i<m;i++)
      {
  		
          if(sol[i]>mx)
          {
              mx = sol[i];
              indx=i;
          }
          //cout<<V[i][t-1]<<" "<<mx<<" "<<indx<<endl;
      }
      
      int q;
      q = result[t-1] = indx;
       //cout<<q<<endl;
       
      for(int i=t-2;i>=0;i--)
      {
          q = result[i] =  argmnt[q][i+1];
      }
      //print result

      printf("The resultant hidden sequence is : \n");
      for (int i = 0; i < t; ++i)
      {
          printf("%d ",result[i]);
      } 
      printf("\n" );
  }

/* mpi is terminated. */ 
      MPI_Finalize();
  /* end program dot3mpi */
      return 0;
}
