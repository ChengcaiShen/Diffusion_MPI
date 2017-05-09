/*==============================================================================
 Name:
  speed_test_diffequ
 Purpose:
  Test the mpi comunication speed for solving diffusion equations.
 Create:
  By Chengcai
 2017-05-04
  Update:
 2017-05-08
 
 =============================================================================*/
#include "mpi.h"
#include <stdio.h>
#include <math.h>

/*------------------------------------------------------------------------------
  global variables
  ----------------------------------------------------------------------------*/
int nnode = 4; // The number of parallel nodes in each direction
int nx = 512, ny = 512;
double xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax=1.0;


/*------------------------------------------------------------------------------
  Main
  ----------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  /* Define mpich relavited */
  int myid, nprocs, nameln;
  char name[MPI_MAX_PROCESSOR_NAME];
  char message[20];
  
  int id_sta, id_end; // The index of nodes in my group
  int id_recv, id_source;
  int myid_relat;
  int tag_x;
  MPI_Status status;
  
  /* Input and Output */
  FILE *fstream;
  char file_output_mesh[255];
  int total_cycle=500;
  double time_start, time_end;
  
  /* Define local variables */
  int i, j, k, ii, jj;
  int itime;
  double dx, dy;
  double iniw = 0.2;
  double time;
  double recv_buf[ny/nnode][nx], send_buf[ny/nnode][nx];
  double te[ny][nx], x[nx], y[ny]; // standard local var
  double whole_arr_x[ny/nnode][nx*nnode], whole_arr_y[ny*nnode][nx/nnode];
  
  /* (0) Initilize MPICH */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  MPI_Get_processor_name(name, &nameln);
  //printf("Hello! Myid is %d. My hostname is \"%s\"\n", myid, name);
  
  /* (1) Initialize variables */
  // Define x-axis, and y-axis
  double x_sta, y_sta;
  dx = (xmax - xmin)/(nnode*nx-1);
  dy = (ymax - ymin)/(nnode*ny-1);
  x_sta = (double)(myid%nnode)*(dx*nx) + xmin;
  y_sta = (double)(myid%nnode)*(dy*nx) + ymin;
  
  for (i=0; i<=nx-1; i++) {
    x[i] = i*dx + x_sta;
  }
  for (j=0; j<=ny-1; j++) {
    y[j] = j*dy + y_sta;
  }
  //printf("myid=%d, x[0]=%f, x[nx-1]=%f, y[0]=%f, y[ny-1]=%f\n",
  //       myid, x[0], x[nx-1], y[0], y[ny-1]);
  
  // Define the Te
  for (i=0; i<=nx-1; i++) {
    for (j=0; j<=ny-1; j++) {
      //te[j][i] = 10.0*exp(-(x[i]/iniw)*(x[i]/iniw)-(y[j]/iniw)*(y[j]/iniw));
      te[j][i] = myid;
    }
  }
  
  /* (2) Prepare the local buffer in x-direction */
  id_end = (myid/nnode)*nnode + nnode - 1;
  id_sta = id_end - (nnode - 1);
  myid_relat = myid - id_sta;
  
  /*--------------------------------------------------------------------------*/
  /* Enter the time-loop*/
  /*--------------------------------------------------------------------------*/
  time_start = MPI_Wtime();
  for (itime=0; itime <= total_cycle; itime++) {
    if (myid == 0) {
      printf("i_cycle = %d/%d\n", itime, total_cycle);
    }
    
    /* (2.0) MPICH send and recv */
    for (jj=0; jj<=nnode-1; jj++) { // jj is the index of piece
      
      // define send buffer.
      for (j=0; j<=ny/nnode-1; j++) {
        for (i=0; i<=nx-1; i++) {
          send_buf[j][i] = te[j+jj*(ny/nnode)][i];
        }
      }
      
      // send
      if (myid_relat != jj) {
        tag_x = myid;
        id_recv = jj + id_sta;
        MPI_Send(&send_buf, nx*ny/nnode, MPI_DOUBLE, id_recv, tag_x, MPI_COMM_WORLD);
        //printf("When jj=%d, id_%d is sending piece %d to id_%d, tag=%d\n", jj, myid, jj, id_recv, tag_x);
      } else {
        
        // receive all nodes in the same row
        for (ii = 0; ii <= nnode-1; ii++) {
          id_source = ii + id_sta;
          tag_x = ii + id_sta;
          // Put the received piece into the whole_arr_x if they are from others
          if (myid_relat != ii) {
            MPI_Recv(&recv_buf, nx*ny/nnode, MPI_DOUBLE, id_source, tag_x, MPI_COMM_WORLD, &status);
            //printf(" ... id_%d received piece %d from id_%d, tag=%d, \n", myid, jj, id_source, tag_x);
            for (j=0; j<=ny/nnode-1; j++) {
              for (i=0; i<=nx-1; i++) {
                whole_arr_x[j][i+(ii-id_sta)*nx] = recv_buf[j][i];
              }
            }
            
          } else {
            /* Put the piece from myself into the whole_arr_x */
            for (j=0; j<=ny/nnode-1; j++) {
              for (i=0; i<=nx-1; i++) {
                whole_arr_x[j][i+myid_relat*nx] = te[j+jj*(ny/nnode)][i];
              }
            }
          }
        } // end receive here
      } // end either send or receive here
      
      /* MPI sync in the current row */
      MPI_Barrier(MPI_COMM_WORLD);
    }// end all pieces (nnode) sending and receiving here
    
    // show MPI results
    //printf("Myid: %d, whole_arr_x:%f,%f,%f,%f\n", myid,
    //       whole_arr_x[0][0],whole_arr_x[0][nx],whole_arr_x[0][2*nx],whole_arr_x[0][3*nx]);
    
    
    /* (2.1) Computation in x-direction in each node */
    // ... solving diffusion equation in x- ...
    for (j=0; j<=ny/nnode-1; j++ ) {
      for (i=0; i<=nx*nnode-1; i++) {
        whole_arr_x[j][i] = whole_arr_x[j][i] + 0.001*myid;
      }
    }
    //printf(">Myid: %d, whole_arr_x:%f,%f,%f,%f\n", myid, whole_arr_x[0][0],whole_arr_x[0][nx],whole_arr_x[0][2*nx],whole_arr_x[0][3*nx]);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    /* (2.2) Return results to origin blocks */
    // enter the jj_th pieces loop
    for (jj=0; jj<=nnode-1; jj++) {
      
      if (jj != myid_relat) {
        
        // Recv the ii_th piece
        tag_x = myid;
        id_source = jj + id_sta;
        MPI_Recv(&recv_buf, nx*ny/nnode, MPI_DOUBLE, id_source, tag_x, MPI_COMM_WORLD, &status);
        //printf(" === id_%d received piece jj.%d from id_%d, tag=%d, \n", myid, jj, id_source, tag_x);
        
        // Return recv_buf to origin array (te2D)
        for (j=0; j<=ny/nnode-1; j++) {
          for (i=0; i<=nx-1; i++) {
            te[j+jj*ny/nnode][i] = recv_buf[j][i];
          }
        }
        
      } else {
        
        // Send the ii_th of whole_arr_x to other nodes
        for (ii=0; ii<=nnode-1; ii++) {
          
          // Define the send buffer
          for (j=0; j<=ny/nnode-1; j++) {
            for (i=0; i<=nx-1; i++) {
              send_buf[j][i] = whole_arr_x[j][i+ii*nx];
            }
          }
          
          // if ii ne jj then sending (the current ii_th piece is not belong myid)
          // else return the ii_th to the origin array (e.g., te2D)
          tag_x = ii + id_sta;
          id_recv = ii + id_sta;
          if (ii != jj) {
            MPI_Send(&send_buf, nx*ny/nnode, MPI_DOUBLE, id_recv, tag_x, MPI_COMM_WORLD);
            //printf("= id_%d is sending the piece ii.%d to id_%d(recv), tag=%d\n", myid, ii, id_recv, tag_x);
          } else {
            for (j=0; j<=ny/nnode-1; j++) {
              for (i=0; i<=nx-1; i++) {
                te[j+jj*ny/nnode][i] = send_buf[j][i];
              }
            }
          }
          
        } // end ii_th loop in here
        
      } // end if Send (or Recv) in here
      //MPI_Barrier(MPI_COMM_WORLD);
      
    } // end jj piece loop in here
    
    
  } // end time-loop inhere
  
  // Get the time_end and speed
  time_end = MPI_Wtime();
  if (myid == 0) {
    printf("Cost = %f(s); Cost/step = %f(s)\n", time_end-time_start, (time_end-time_start)/(total_cycle+1));
  }
  
  /* (3) Write results into files */
  time = 0.0;
  sprintf(file_output_mesh,"id%04d_t%06.2f.dat", myid, time);
  //printf("Write resultes to file: %s\n", file_output_mesh);
  fstream = fopen(file_output_mesh,"w");
  fprintf(fstream, "%d\n", nx);
  fprintf(fstream, "%d\n", ny);
  fprintf(fstream, "%f\n", time);
  for(i = 0; i <= nx-1; i++) {
    fprintf(fstream, "%f\n", x[i]);
  }
  for(j = 0; j <= ny-1; j++) {
    fprintf(fstream, "%f\n", y[j]);
  }
  for(i = 0; i <= nx-1; i++) {
    for(j = 0; j <= ny-1; j++) {
      fprintf(fstream, "%f\n", te[j][i]);
    }
  }
  fclose(fstream);
  
  /* Normal stop. */
  MPI_Finalize();
  return(0);
}

