/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_IRREGULAR_H
#define SPK_IRREGULAR_H

#include "pointers.h"

namespace SPPARKS_NS {

class Irregular : protected Pointers {
 public:

  // static variable across all Irregular objects, for qsort callback

  static int *proc_recv_copy;

  Irregular(class SPPARKS *);
  ~Irregular();
  int create_data(int, int *, int sortflag = 0);
  void exchange_data(char *, int, char *);
  void destroy_data();
  bigint memory_usage();

 private:
  int me,nprocs;

  int maxsend,maxrecv;              // size of buf send/recv in # of doubles
  double *buf_send,*buf_recv;       // bufs used in migrate_atoms
  int maxdbuf;                      // size of double buf in bytes
  double *dbuf;                     // double buf for largest single atom send
  int maxbuf;                       // size of char buf in bytes
  char *buf;                        // char buf for largest single data send

  int *mproclist,*msizes;           // persistent vectors in migrate_atoms
  int maxlocal;                     // allocated size of mproclist and msizes

  int *work1,*work2;                // work vectors

  // plan params for irregular communication of atoms or datums
  // no params refer to atoms/data copied to self

  int nsend_proc;            // # of messages to send
  int nrecv_proc;            // # of messages to recv
  int sendmax_proc;          // # of doubles/datums in largest send message
  int *proc_send;            // list of procs to send to
  int *num_send;             // # of atoms/datums to send to each proc
  int *index_send;           // list of which atoms/datums to send to each proc
  int *proc_recv;            // list of procs to recv from
  MPI_Request *request;      // MPI requests for posted recvs
  MPI_Status *status;        // MPI statuses for WaitAll

  // extra plan params plan for irregular communication of atoms
  // no params refer to atoms copied to self

  int *length_send;          // # of doubles to send to each proc
  int *length_recv;          // # of doubles to recv from each proc
  int *offset_send;          // where each atom starts in send buffer

  // extra plan params plan for irregular communication of datums
  // 2 self params refer to data copied to self

  int *num_recv;             // # of datums to recv from each proc
  int num_self;              // # of datums to copy to self
  int *index_self;           // list of which datums to copy to self

  // private methods

  int create_atom(int, int *, int *, int);
  void exchange_atom(double *, int *, double *);
  void destroy_atom();

  int binary(double, int, double *);

  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
};

}

#endif

/* ERROR/WARNING messages:

*/
