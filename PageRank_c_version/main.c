#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>

// This block enables compilation of the code with and without LIKWID in place
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif




void csr_spmv(int m, int const * r, int const * j, float const * a, float const * x, float * y);
void coo_spmv(int N, int const * i, int const * j, float const * a, float const * x, float * y);




int main(int argc, char *argv[])
{
      LIKWID_MARKER_INIT;
      /*************************** TIME, VARIABLES ***************************/
      
      // Keep track of the execution time
      clock_t begin, end;
      double time_spent;
      begin = clock();

      /******************* OPEN FILE + NUM OF NODES/EDGES ********************/

      // Open the data set
      char *filename = argv[1];
      
      FILE *fp;
      if((fp = fopen(filename,"r")) == NULL)
      {
        fprintf(stderr,"[Error] Cannot open the file");
        exit(1);
      }

      // Read the data set and get the number of nodes (n) and edges (e)
      int n, e;
      char ch;
      char str[100];
      ch = getc(fp);
      while(ch == '#') 
      {
        fgets(str,100-1,fp);
        sscanf (str,"%*s %d %*s %d", &n, &e); //number of nodes
        ch = getc(fp);
      }
      ungetc(ch, fp);
      
      // DEBUG: Print the number of nodes and edges, skip everything else
      printf("\nGraph data:\n\n  Nodes: %d, Edges: %d \n\n", n, e);
      
      
      /************************* CSR STRUCTURES *****************************/
        
      /* Compressed sparse row format: 
        - Val vector: contains 1.0 if an edge exists in a certain row
        - Col_ind vector: contains the column index of the corresponding value in 'val'
        - Row_ptr vector: points to the start of each row in 'col_ind'
      */

      
      float *val = calloc(e, sizeof(float));
      int *col_ind = calloc(e, sizeof(int));
      int *row_ind = calloc(e, sizeof(int));
      int *row_ptr = calloc(n+1, sizeof(int));


      
      // The first row always starts at position 0
      row_ptr[0] = 0;
      
      int fromnode, tonode;
      int cur_row = 0;
      int i = 0;
      int j = 0;
      // Elements for row
      int elrow = 0;
      // Cumulative numbers of elements
      int curel = 0;
      
      
      
      while(!feof(fp))
      {
        
          fscanf(fp,"%d%d",&fromnode,&tonode);
          if (fromnode > cur_row) 
          { 
            // change the row
            curel = curel + elrow;
            for (int k = cur_row + 1; k <= fromnode; k++) 
            {
              row_ptr[k] = curel;
            }
            elrow = 0;
            cur_row = fromnode;
          }
          val[i] = 1.0;
          col_ind[i] = tonode;
          row_ind[i] = cur_row;
          elrow++;
          i++;
      }

      row_ptr[cur_row+1] = curel + elrow - 1;


      // Fix the stochastization
      int *out_link = (int*) malloc(n*sizeof(int));
      for(i=0; i<n; i++)
      {
        out_link[i] =0;
      }

      int rowel = 0;
      for(i=0; i<n; i++)
      {
            if (row_ptr[i+1] != 0) 
            {
              rowel = row_ptr[i+1] - row_ptr[i];
              out_link[i] = rowel;
            }
      }
        
      int curcol = 0;
      for(i=0; i<n; i++)
      {
        rowel = row_ptr[i+1] - row_ptr[i];
        for (j=0; j<rowel; j++) 
        {
          val[curcol] = val[curcol] / out_link[i];
          curcol++;
        }
      }


    
      /******************* INITIALIZATION OF P, DAMPING FACTOR ************************/

      // Set the damping factor 'd'
      float d = 0.85;
      
      // Initialize p[] vector
      float *p = (float*) malloc(n*sizeof(float));
      for(i=0; i<n; i++)
      {
        p[i] = 1.0/n;
      }
      
      /*************************** PageRank LOOP  **************************/

      // Set the looping condition and the number of iterations 'k'
      int looping = 1;
      int k = 0;
      
      // Initialize new p vector
      float *p_new = (float*) malloc(n*sizeof(float));
      
      double itime, ftime, exec_time;
      itime = omp_get_wtime();
      while (looping){

      // Initialize p_new as a vector of n 0.0 cells
      for(i=0; i<n; i++)
      {
        p_new[i] = 0.0;
      }
        
      int rowel = 0;
      int curcol = 0;
      
      //csr_spmvcsr_spmv(n, row_ptr, col_ind, val, p, p_new);
      coo_spmv(e, row_ind, col_ind, val, p, p_new);

        // Page rank modified algorithm 
       /*    LIKWID_MARKER_START("pageRank");
        for(i=0; i<n; i++)
        {
          rowel = row_ptr[i+1] - row_ptr[i];
          for (j=0; j<rowel; j++) 
          {
            p_new[col_ind[curcol]] = p_new[col_ind[curcol]] + val[curcol] * p[i];
            curcol++;
          }
        }
        LIKWID_MARKER_STOP("pageRank");*/
        // Adjustment to manage dangling elements 
        for(i=0; i<n; i++)
        {
          p_new[i] = d * p_new[i] + (1.0 - d) / n;
        }
          
        // TERMINATION: check if we have to stop
        float error = 0.0;
        for(i=0; i<n; i++) 
        {
          error =  error + fabs(p_new[i] - p[i]);
        }
        //if two consecutive instances of pagerank vector are almost identical, stop
        if (error < 0.000001)
        {
          looping = 0;
        }
        
        // Update p[]
        for (i=0; i<n;i++)
        {
            p[i] = p_new[i];
        }
        
        // Increase the number of iterations
        k = k + 1;
    }
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("\n\nTime of parallel is : %f", exec_time);
      
    /*************************** CONCLUSIONS *******************************/

    // Stop the timer and compute the time spent
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

      
    // Print results
    printf ("\nNumber of iteration to converge: %d \n\n", k);
    printf ("\nsize of node: %d \n\n", n); 
    /*printf ("Final Pagerank values:\n\n[");
    for (i=0; i<n; i++){
      printf("%f ", p[i]);
      if(i!=(n-1)){ printf(", "); }
    } */
    printf("]\n\nTime spent: %f seconds.\n", time_spent);
    LIKWID_MARKER_CLOSE;
    return 0;
}



void csr_spmv(int m, int const * r, int const * j, float const * a, float const * x, float * y) 
{
  #pragma omp parallel
  {
    LIKWID_MARKER_START("pageRank");
    #pragma omp for schedule(static)
    for (int i = 0; i < m; i++) 
    {
      float z = 0.0;
      for (int k = r[i]; k < r[i+1]; k++)
      {
        z += a[k] * x[j[k]]; 
      }
      y[i] += z;
    }
    LIKWID_MARKER_STOP("pageRank");
  }
}


void coo_spmv(int N, int const * i, int const * j, float const * a, float const * x, float * y) 
{
  #pragma omp parallel
  {
    LIKWID_MARKER_START("pageRank");
    #pragma omp for schedule(static)
    for ( int k = 0; k < N; k++) 
    {
      #pragma omp atomic
      y[i[k]] += a[k] * x[j[k]];
    }
    LIKWID_MARKER_STOP("pageRank");
  }
}

