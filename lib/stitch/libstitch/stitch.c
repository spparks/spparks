/* Copyright 2019 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525
 * with NTESS, the U.S. Government retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * For more information, contact Jay Lofstead (gflofst@sandeia.gov) or
 * John Mitchell (jamitch@sandia.gov) for more information.
 */ 
// use this define to get strdup defined and some other functions
#define _SVID_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <inttypes.h>
#include <errno.h>
#ifdef STITCH_TIMING
// only needed if timing is collected
#include <sys/time.h>
#endif

#ifdef STITCH_PARALLEL
#include "mpi.h"
#endif

#include "stitch.h"
#include "sqlite3.h"

enum
{
     STITCH_MPI_TAG_WRITE_BLOCK = 1000
    ,STITCH_MPI_TAG_READ_BLOCK = 1001
    ,STITCH_MPI_TAG_SET_PARAMETERS_BLOCK = 1002
};

typedef struct StitchFile_struct
{
    sqlite3 * db;
#ifdef STITCH_PARALLEL
    MPI_Comm comm;
    int size;
#endif
    int rank;
} StitchFile;


static int callback (void * NotUsed, int argc, char ** argv, char ** ColName)
{
    int i;
    for (i = 0; i < argc; i++)
    {
        printf ("%s = %s\n", ColName [i], argv [i] ? argv [i] : "NULL");
    }
    printf ("\n");
    return 0;
}

#ifdef STITCH_PARALLEL
int stitch_open (StitchFile ** file, MPI_Comm comm, const char * filename)
#else
int stitch_open (StitchFile ** file, const char * filename)
#endif
{
    int rc = 0;
    char * ErrMsg = 0;

    *file = (StitchFile *) malloc (sizeof (StitchFile));
    memset (*file, 0, sizeof (StitchFile)); // sets rank to 0 by default

#ifdef STITCH_PARALLEL
    rc = MPI_Comm_dup (comm, &(*file)->comm);
    assert (rc == MPI_SUCCESS);

    int rank = 0;
    MPI_Comm_rank ((*file)->comm, &rank);
    (*file)->rank = rank;
    MPI_Comm_size ((*file)->comm, &((*file)->size));

    // only rank 0 does the open and db create
    if (rank == 0)
    {
#endif
        // check to see if file exists already
        struct stat f;
        bool file_exists = false;
        rc = stat (filename, &f);
        if (rc == -1)
        {
            if (errno != ENOENT)
            {
                fprintf (stderr, "error open/creating file %s: %d\n", filename, errno);
                goto cleanup;
            }
            else
            {
                //printf ("file %s does not exist yet: %d\n", filename, errno);
            }
        }
        else
        {
            //printf ("file %s exists\n", filename);
            file_exists = true;
        }

        //rc = sqlite3_open_v2 (filename, &(*file)->db, SQLITE_OPEN_SHAREDCACHE, 0);
        rc = sqlite3_open_v2 (filename, &(*file)->db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, 0);

        if (rc != SQLITE_OK)
        {
            fprintf (stderr, "Can't open stitch file: %s\n", sqlite3_errmsg ((*file)->db));
            sqlite3_close ((*file)->db);
            goto cleanup;
        }

        if (!file_exists)
        {
            rc = sqlite3_exec ((*file)->db, "create table globals (absolute_tolerance real, relative_tolerance real, no_value_present int, first_time real, last_time real)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            sqlite3_stmt * stmt_index = NULL;
            const char * tail_index = NULL;

            rc = sqlite3_prepare_v2 ((*file)->db, "insert into globals (absolute_tolerance, relative_tolerance, no_value_present) values (?, ?, ?)", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg ((*file)->db));
                sqlite3_close ((*file)->db);
                goto cleanup;
            }
            rc = sqlite3_bind_double (stmt_index, 1, STITCH_DEFAULT_ABSOLUTE_TOLERANCE); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_double (stmt_index, 2, STITCH_DEFAULT_RELATIVE_TOLERANCE); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_int64 (stmt_index, 3, STITCH_NO_VALUE); assert (rc == SQLITE_OK);
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg ((*file)->db));
                sqlite3_close ((*file)->db);
                goto cleanup;
            }
            assert (rc == SQLITE_OK || rc == SQLITE_DONE);
            rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK);

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            rc = sqlite3_exec ((*file)->db, "create table times (timestamp integer primary key, time real not null)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index times_time on times (time)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create unique index times_timestamp on times (timestamp)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            rc = sqlite3_exec ((*file)->db, "create table fields (field_id integer primary key, field text not null collate nocase, type int not null, length int not null, no_value_int int, no_value_real real)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create unique index fields_field_id on fields (field_id)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            rc = sqlite3_exec ((*file)->db, "create table blocks (field_id int not null, timestamp int not null, x_min int not null, y_min int not null, z_min int not null, x_max int not null, y_max int not null, z_max int not null, state blob)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index field_id_blocks on blocks (field_id)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_timestamp on blocks (timestamp)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_x_min on blocks (x_min)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_y_min on blocks (y_min)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_z_min on blocks (z_min)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_x_max on blocks (x_max)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_y_max on blocks (y_max)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }

            rc = sqlite3_exec ((*file)->db, "create index blocks_z_max on blocks (z_max)", callback, 0, &ErrMsg);
            if (rc != SQLITE_OK)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, ErrMsg);
                sqlite3_free (ErrMsg);
                sqlite3_close ((*file)->db);
                goto cleanup;
            }
        }
#ifdef STITCH_PARALLEL
        MPI_Barrier ((*file)->comm);
    }
    else
    {
        // be sure to wait until the file is created properly before everyone else opens it
        MPI_Barrier ((*file)->comm);
        //rc = sqlite3_open_v2 (filename, &(*file)->db, SQLITE_OPEN_SHAREDCACHE, 0);
        rc = sqlite3_open_v2 (filename, &(*file)->db, SQLITE_OPEN_READWRITE, 0);

        if (rc != SQLITE_OK)
        {
            fprintf (stderr, "%d Can't open stitch file: %s\n", (*file)->rank, sqlite3_errmsg ((*file)->db));
            sqlite3_close ((*file)->db);
            goto cleanup;
        }
    }
#endif

    // this should reduce pressure on /tmp and improve performance. If memory is tight, this might
    // be a problem.
    rc = sqlite3_exec ((*file)->db, "PRAGMA temp_store = MEMORY", NULL, NULL, &ErrMsg);
    if (rc != SQLITE_OK)
    {
        fprintf (stderr, "%d Can't open stitch file: %s\n", (*file)->rank, sqlite3_errmsg ((*file)->db));
        sqlite3_close ((*file)->db);
        goto cleanup;
    }

cleanup:
    return rc;
}

int stitch_close (StitchFile ** file)
{
    sqlite3_close ((*file)->db);
#ifdef STITCH_PARALLEL
    MPI_Comm_free (&(*file)->comm);
#endif
    memset (*file, 0, sizeof (StitchFile));
    free (*file);
    *file = 0;

    return 0;
}

int stitch_set_parameters (const StitchFile * file, double absolute_tolerance, double relative_tolerance, int64_t default_no_value_present)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
#ifdef STITCH_PARALLEL
    MPI_Status status;

    if (file->rank == 0)
    {
#endif
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select absolute_tolerance from globals", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_DONE && rc != SQLITE_ROW)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

        if (rc == SQLITE_DONE)
        {
            rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK);
            do
            {
                rc = sqlite3_prepare_v2 (file->db, "insert into globals (absolute_tolerance, relative_tolerance, no_value_present) values (?, ?, ?)", -1, &stmt_index, &tail_index);
                if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

            rc = sqlite3_bind_double (stmt_index, 1, absolute_tolerance); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_double (stmt_index, 2, relative_tolerance); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_int64 (stmt_index, 3, default_no_value_present); assert (rc == SQLITE_OK);

            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        }
        else
        {
            if (rc == SQLITE_ROW)
            {
                do
                {
                    rc = sqlite3_finalize (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
                do
                {
                    rc = sqlite3_prepare_v2 (file->db, "update globals set absolute_tolerance = ?, relative_tolerance = ?, no_value_present = ?", -1, &stmt_index, &tail_index);
                    if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

                rc = sqlite3_bind_double (stmt_index, 1, absolute_tolerance); assert (rc == SQLITE_OK);
                rc = sqlite3_bind_double (stmt_index, 2, relative_tolerance); assert (rc == SQLITE_OK);
                rc = sqlite3_bind_int64 (stmt_index, 3, default_no_value_present); assert (rc == SQLITE_OK);

                do
                {
                    rc = sqlite3_step (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
            }
            else
            {
                fprintf (stderr, "Line: %d what? rc = %d\n", __LINE__, rc);
                assert (0);
            }
        }

        rc = sqlite3_finalize (stmt_index);
#ifdef STITCH_PARALLEL
//        printf ("[%d] send to %d\n", file->rank, file->rank+1);
        MPI_Send (&rc, 1, MPI_INT, file->rank + 1, STITCH_MPI_TAG_SET_PARAMETERS_BLOCK, file->comm);
//        printf ("[%d] send done\n");
    }
    else
    {
//        printf ("[%d] recv from %d\n", file->rank, file->rank-1);
        MPI_Recv (&rc, 1, MPI_INT, file->rank - 1, STITCH_MPI_TAG_SET_PARAMETERS_BLOCK, file->comm, &status);
        if (file->rank < file->size - 1)
            MPI_Send (&rc, 1, MPI_INT, file->rank + 1, STITCH_MPI_TAG_SET_PARAMETERS_BLOCK, file->comm);
//        printf ("[%d] recv done\n", file->rank);
    }
#endif
//    int rowid = sqlite3_last_insert_rowid (file->db);
//printf ("rowid for chunk: %d\n", rowid);

//printf ("insert chunk completed\n");
#ifdef STITCH_PARALLEL
//    MPI_Barrier (file->comm); // make sure all is set before anyone continues
#endif

cleanup:

     return rc;
}

int stitch_get_parameters (const StitchFile * file, double * absolute_tolerance, double * relative_tolerance, int64_t * default_no_value_present, double * first_time, double * last_time)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

    do
    {
        rc = sqlite3_prepare_v2 (file->db, "select absolute_tolerance, relative_tolerance, no_value_present, first_time, last_time from globals", -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf (stderr, "%d get_parameters, prepare Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

    do
    {
        rc = sqlite3_step (stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf (stderr, "%d get_parameters, prepare Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

    while (rc == SQLITE_ROW)
    {
        *absolute_tolerance = sqlite3_column_double (stmt_index, 0);
        *relative_tolerance = sqlite3_column_double (stmt_index, 1);
        *default_no_value_present = sqlite3_column_int64 (stmt_index, 2);
        *first_time = sqlite3_column_double (stmt_index, 3);
        *last_time = sqlite3_column_double (stmt_index, 4);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_DONE && rc != SQLITE_ROW)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    }

cleanup:
    rc = sqlite3_finalize (stmt_index);
    if (rc != SQLITE_OK)
    {
        fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
        sqlite3_close (file->db);
    }

    return 0;
}

static int stitch_create_timestamp (const StitchFile * file, double * time, int64_t * timestamp, int32_t * is_new_time)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
    double real_time = *time;

#ifdef STITCH_PARALLEL
    if (file->rank == 0)
    {
#endif
        do
        {
            // check to see if the globals was initialized. If not, default to 0.0 tolerance
            rc = sqlite3_prepare_v2 (file->db, "select * from globals", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        rc = sqlite3_step (stmt_index);

        // if DONE, then no initialization done
        if (rc == SQLITE_DONE)
        {
            sqlite3_finalize (stmt_index);
            do
            {
                rc = sqlite3_prepare_v2 (file->db, "insert into globals (absolute_tolerance, relative_tolerance) values (?, ?)", -1, &stmt_index, &tail_index);
                if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
            rc = sqlite3_bind_double (stmt_index, 1, STITCH_DEFAULT_ABSOLUTE_TOLERANCE); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_double (stmt_index, 2, STITCH_DEFAULT_RELATIVE_TOLERANCE); assert (rc == SQLITE_OK);
            do
            {
                rc = sqlite3_step (stmt_index); assert (rc == SQLITE_OK);
                if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
        }
        rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK);

        //printf ("create_timestamp: time: %g\n", real_time);
        // check to see if the time exists already. If so, use that timestamp.
        // otherwise, create a new timestamp.
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select timestamp, time from times, globals where abs(times.time - ?) < (globals.absolute_tolerance + (? * globals.relative_tolerance))", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
        rc = sqlite3_bind_double (stmt_index, 1, real_time); assert (rc == SQLITE_OK);
        rc = sqlite3_bind_double (stmt_index, 2, real_time); assert (rc == SQLITE_OK);

        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
        if (rc == SQLITE_ROW)
        {
            // get the timestamp and use that
            *timestamp = sqlite3_column_int (stmt_index, 0);
            *time = sqlite3_column_double (stmt_index, 1);
            //printf ("create_timestamp: got row timestamp: %" PRId64 " time: %g\n", *timestamp, *time);
            // note that we got an existing time
            *is_new_time = false;
            rc = sqlite3_finalize (stmt_index);
        }
        else
        {
            // close out the query before starting the insert
            sqlite3_finalize (stmt_index);

            if (rc == SQLITE_DONE)
            {
                // no time, so create a new one
                // note that we created a new time
                *is_new_time = true;
                do
                {
                    rc = sqlite3_prepare_v2 (file->db, "insert into times (time) values (?)", -1, &stmt_index, &tail_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
                if (rc != SQLITE_OK)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }

                rc = sqlite3_bind_double (stmt_index, 1, *time); assert (rc == SQLITE_OK);

                do
                {
                    rc = sqlite3_step (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_ROW  && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_LOCKED)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

                rc = sqlite3_finalize (stmt_index);
                *timestamp = sqlite3_last_insert_rowid (file->db);
                //printf ("create_timestamp: no row. adding %g timestamp: %" PRId64 "\n", *time, *timestamp);

                do
                {
                    rc = sqlite3_prepare_v2 (file->db, "select first_time, last_time from globals", -1, &stmt_index, &tail_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

                do
                {
                    rc = sqlite3_step (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

                int is_null_first = sqlite3_column_type (stmt_index, 0);

                double first_time = *time;
                double last_time = *time;
                if (is_null_first != SQLITE_NULL)
                {
                    first_time = sqlite3_column_double (stmt_index, 0);
                    last_time = sqlite3_column_double (stmt_index, 1);

                    if (*time < first_time)
                        first_time = *time;
                    if (*time > last_time)
                        last_time = *time;
                }
                rc = sqlite3_finalize (stmt_index);

                do
                {
                    rc = sqlite3_prepare_v2 (file->db, "update globals set first_time = ?, last_time = ?", -1, &stmt_index, &tail_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
                rc = sqlite3_bind_double (stmt_index, 1, first_time);
                rc = sqlite3_bind_double (stmt_index, 2, last_time);
                do
                {
                    rc = sqlite3_step (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

                rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK);
            }
            else
            {
                fprintf (stderr, "Line: %d rc: %d\n", __LINE__, rc);
                assert (0);
                // what's going on here?
            }
        }
#ifdef STITCH_PARALLEL
    }
    //printf ("1 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (timestamp, 1, MPI_INT64_T, 0, file->comm);
#endif

cleanup:
    return rc;
}

static int stitch_get_timestamp (const StitchFile * file, double * time, int64_t * timestamp, int32_t * is_new_time)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

#ifdef STITCH_PARALLEL
    if (file->rank == 0)
    {
#endif
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select timestamp, time from times, globals where abs(times.time - ?) < (globals.absolute_tolerance + (? * globals.relative_tolerance))", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_LOCKED && rc != SQLITE_ROW)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

        rc = sqlite3_bind_double (stmt_index, 1, *time); assert (rc == SQLITE_OK);
        rc = sqlite3_bind_double (stmt_index, 2, *time); assert (rc == SQLITE_OK);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

        if (rc == SQLITE_ROW)
        {
            *timestamp = sqlite3_column_int (stmt_index, 0);
            *time = sqlite3_column_double (stmt_index, 1);
            *is_new_time = false;
        }
        else
        {
            if (rc == SQLITE_DONE)
            {
                *timestamp = -1;
                *time = 0;
                *is_new_time = true;
            }
            else
            {
                *timestamp = -1;
                *time = 0;
                *is_new_time = true;
                // what's going on here?
                fprintf (stderr, "Line: %d rc: %d\n", __LINE__, rc);
                rc = sqlite3_finalize (stmt_index);
                sqlite3_close (file->db);
                assert (0);
            }
        }
        rc = sqlite3_finalize (stmt_index);
#ifdef STITCH_PARALLEL
    }
    //printf ("2 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (timestamp, 1, MPI_INT64_T, 0, file->comm);
    //printf ("3 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (time, 1, MPI_DOUBLE, 0, file->comm);
    //printf ("4 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (is_new_time, 1, MPI_INT32_T, 0, file->comm);
#endif

    //printf ("rank: %d timestamp: %" PRId64 "\n", file->rank, *timestamp);

cleanup:

    return rc;
}

int stitch_get_times_count (const StitchFile * file, int64_t * count)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

#ifdef STITCH_PARALLEL
    if (file->rank == 0)
    {
#endif
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select count (*) from times", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        rc = sqlite3_step (stmt_index); //assert (rc == SQLITE_ROW || rc == SQLITE_DONE);

        while (rc == SQLITE_ROW)
        {
            *count = sqlite3_column_int (stmt_index, 0);
            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
        }

cleanup:
        sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK || rc == SQLITE_DONE);
#ifdef STITCH_PARALLEL
    }

    //printf ("5 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (count, 1, MPI_INT64_T, 0, file->comm);
#endif

    return rc;
}

int stitch_get_times (const StitchFile * file, int64_t * count, double * times)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
    int64_t buffer_entries = *count;
    *count = 0; // how many we retrieved. On input, it is how big the buffer is

#ifdef STITCH_PARALLEL
    if (file->rank == 0)
    {
#endif
        // strictly speaking, ordering by the time is more correct since SQLite might reuse a deleted timestamp entry
        // given our potential output size, this should NEVER be an issue. It is several orders of magnitiude beyond
        // the most extreme case.  Ordering by an integer should be faster than a float
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select time from times order by time", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        assert (stmt_index != NULL);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

        while (rc == SQLITE_ROW && *count < buffer_entries)
        {
//            timestamps [*count] = sqlite3_column_int64 (stmt_index, 0);
            times [*count]  = sqlite3_column_double (stmt_index, 0);
            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
            (*count)++;
        }

cleanup:
        sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK || rc == SQLITE_DONE);
#ifdef STITCH_PARALLEL
    }

    //printf ("6 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (count, 1, MPI_INT64_T, 0, file->comm);
    //MPI_Bcast (timestamps, *count, MPI_INT64_T, 0, file->comm);
    //printf ("7 MPI_Bcast %d\n", file->rank);
    MPI_Bcast (times, *count, MPI_DOUBLE, 0, file->comm);
#endif

    return rc;
}

// retrieve the number of fields
int stitch_get_fields_count (const StitchFile * file, int64_t * count)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

    if (file->rank == 0)
    {
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select count (*) from fields", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        while (rc == SQLITE_ROW)
        {
            *count = sqlite3_column_int (stmt_index, 0);
            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        }

cleanup:
        sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK || rc == SQLITE_DONE);
    }

#ifdef STITCH_PARALLEL
    //printf ("7a MPI_Bcast %d\n", file->rank);
    MPI_Bcast (count, 1, MPI_INT64_T, 0, file->comm);
#endif

    return rc;
}

// count is an in/out parameter. For IN, it is the size of the buffers. For OUT, it is the number of entries.
// the labels are allocated with malloc and must be freed by the caller.
int stitch_get_fields (const StitchFile * file, int64_t * count, int64_t * field_ids, char ** labels, enum STITCH_TYPES * types, int32_t * lengths, union StitchTypesUnion * no_value_present_values)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
    int64_t buffer_entries = *count;
    *count = 0; // how many we retrieved. On input, it is how big the buffer is

// should do this as a single proc and broadcast, but broadcasting an array of strings is a pain. Just hammer the file for now.
//    if (file->rank == 0)
    {
        // strictly speaking, ordering by the time is more correct since SQLite might reuse a deleted timestamp entry
        // given our potential output size, this should NEVER be an issue. It is several orders of magnitiude beyond
        // the most extreme case.  Ordering by an integer should be faster than a float
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select field_id, field, type, length, no_value_int, no_value_real from fields order by field_id", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        assert (stmt_index != NULL);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        while (rc == SQLITE_ROW && *count < buffer_entries)
        {
            field_ids [*count]  = sqlite3_column_int (stmt_index, 0);
            labels [*count] = strdup ((char *) sqlite3_column_text (stmt_index, 1));
            types [*count]  = sqlite3_column_int (stmt_index, 2);
            lengths [*count]  = sqlite3_column_int (stmt_index, 3);
            switch (types [*count])
            {
                case STITCH_INT32:
                    no_value_present_values [*count].i32 = sqlite3_column_int (stmt_index, 4);
                    break;

                case STITCH_INT64:
                    no_value_present_values [*count].i64 = sqlite3_column_int64 (stmt_index, 4);
                    break;

                case STITCH_FLOAT64:
                    no_value_present_values [*count].f64 = sqlite3_column_double (stmt_index, 5);
                    break;

                case STITCH_NO_TYPE:
                default:
                    fprintf (stderr, "Line: %d. unknown type: %d.\n", __LINE__, types [*count]);
            }
            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
            (*count)++;
        }

cleanup:
        rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK || rc == SQLITE_DONE);
    }

// should do this as a broadcast, but broadcasting an array of strings is a pain. Just doing everyone hammer the server at the moment
    //printf ("6a MPI_Bcast %d\n", file->rank);
//    MPI_Bcast (count, 1, MPI_INT64_T, 0, file->comm);
    //MPI_Bcast (timestamps, *count, MPI_INT64_T, 0, file->comm);
    //printf ("7b MPI_Bcast %d\n", file->rank);
//    MPI_Bcast (times, *count, MPI_DOUBLE, 0, file->comm);

    return rc;
}

// the Python interface needs to determine which of the type-safe C functions to call based on the type of the field being written.
int stitch_get_field_type (const StitchFile * file, int64_t field_id, enum STITCH_TYPES * type)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

    //if (file->rank == 0)
    {
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select type from fields where field_id = ?", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        rc = sqlite3_bind_int (stmt_index, 1, field_id);
        do
        {
            rc = sqlite3_step (stmt_index); //assert (rc == SQLITE_ROW || rc == SQLITE_DONE);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        while (rc == SQLITE_ROW)
        {
            *type = sqlite3_column_int (stmt_index, 0);
            do
            {
                rc = sqlite3_step (stmt_index);
                if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
            } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        }

cleanup:
        rc = sqlite3_finalize (stmt_index); assert (rc == SQLITE_OK || rc == SQLITE_DONE);
    }

    //printf ("5 MPI_Bcast %d\n", file->rank);
//    MPI_Bcast (count, 1, MPI_INT64_T, 0, file->comm);

    return rc;
}

// return value is a ID that corresponds to this field in the file (field_id).
// If the field does not exist yet, create it. Otherwise, reuse it.
int stitch_create_field (const StitchFile * file, const char * label, enum STITCH_TYPES type, union StitchTypesUnion no_value_present, int32_t per_site_length, int64_t * field_id)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

#ifdef STITCH_PARALLEL
    if (file->rank == 0)
#endif
    {
        assert (type == STITCH_INT32 || type == STITCH_INT64 || type == STITCH_FLOAT64);

        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select field_id from fields where field = ?", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        rc = sqlite3_bind_text (stmt_index, 1, strdup (label), strlen (label), free); assert (rc == SQLITE_OK);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        if (rc == SQLITE_ROW)
        {
            //uint32_t db_type = 0;
            //int32_t db_length = 0;
            //const unsigned char * name = NULL;

            *field_id = sqlite3_column_int64 (stmt_index, 0);
            rc = sqlite3_finalize (stmt_index);
        }
        else
        {
            if (rc == SQLITE_DONE)
            {
                //printf ("found nothing\n");
                // close out the query before the insert
                rc = sqlite3_finalize (stmt_index);

                do
                {
                    rc = sqlite3_prepare_v2 (file->db, "insert into fields (field, type, length, no_value_int, no_value_real) values (?, ?, ?, ?, ?)", -1, &stmt_index, &tail_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
                if (rc != SQLITE_OK)
                {
                    fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                    sqlite3_close (file->db);
                    goto cleanup;
                }
                rc = sqlite3_bind_text (stmt_index, 1, strdup (label), strlen (label), free); assert (rc == SQLITE_OK);
                rc = sqlite3_bind_int (stmt_index, 2, type); assert (rc == SQLITE_OK);
                rc = sqlite3_bind_int (stmt_index, 3, per_site_length); assert (rc == SQLITE_OK);
                switch (type)
                {
                    case STITCH_INT32:
                        //printf ("new field no value int32: %d\n", *(int32_t *) no_value_present);
                        rc = sqlite3_bind_int (stmt_index, 4, no_value_present.i32); assert (rc == SQLITE_OK);
                        rc = sqlite3_bind_null (stmt_index, 5); assert (rc == SQLITE_OK);
                        break;

                    case STITCH_INT64:
                        //printf ("new field no value int64: %lld\n", *(int64_t *) no_value_present);
                        rc = sqlite3_bind_int64 (stmt_index, 4, no_value_present.i64); assert (rc == SQLITE_OK);
                        rc = sqlite3_bind_null (stmt_index, 5); assert (rc == SQLITE_OK);
                        break;

                    case STITCH_FLOAT64:
                        //printf ("new field no value int64: %g\n", *(double *) no_value_present);
                        rc = sqlite3_bind_null (stmt_index, 4); assert (rc == SQLITE_OK);
                        rc = sqlite3_bind_double (stmt_index, 5, no_value_present.f64); assert (rc == SQLITE_OK);
                        break;

                    case STITCH_NO_TYPE:
                    default:
                        fprintf (stderr, "Line: %d invalid field type found: %d\n", __LINE__, type);
                        assert (0);
                        break;
                }

                do
                {
                    rc = sqlite3_step (stmt_index);
                    if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
                    {
                        fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                        sqlite3_close (file->db);
                        goto cleanup;
                    }
                } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

                rc = sqlite3_finalize (stmt_index);
                *field_id = sqlite3_last_insert_rowid (file->db);
            }
            else
            {
                *field_id = -1;
                // what's going on here?
                fprintf (stderr, "Line: %d rc: %d\n", __LINE__, rc);
                rc = sqlite3_finalize (stmt_index);
                sqlite3_close (file->db);
                assert (0);
            }
        }
        //printf ("rank: %d new field_id: %" PRId64 "\n", file->rank, *field_id);
    }
#ifdef STITCH_PARALLEL
    MPI_Bcast (field_id, 1, MPI_INT64_T, 0, file->comm);
#endif

cleanup:

    return rc;
}

int stitch_set_field_no_value_present (const StitchFile * file, int64_t field_id, union StitchTypesUnion no_value_present)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
    enum STITCH_TYPES type;

    rc = stitch_get_field_type (file, field_id, &type);

    do
    {
        rc = sqlite3_prepare_v2 (file->db, "update fields set int_no_value = ?, real_no_value = ? where field_id = ?", -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK || rc != SQLITE_BUSY || rc != SQLITE_LOCKED)
        {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
        }
    } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);

    switch (type)
    {
        case STITCH_INT32:
            //printf ("new field no value int32: %d\n", *(int32_t *) no_value_present);
            rc = sqlite3_bind_int (stmt_index, 1, no_value_present.i32); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_null (stmt_index, 2); assert (rc == SQLITE_OK);
            break;

        case STITCH_INT64:
            //printf ("new field no value int64: %" PRId64 "\n", *(int64_t *) no_value_present);
            rc = sqlite3_bind_int64 (stmt_index, 1, no_value_present.i64); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_null (stmt_index, 2); assert (rc == SQLITE_OK);
            break;

        case STITCH_FLOAT64:
            //printf ("new field no value int64: %g\n", *(double *) no_value_present);
            rc = sqlite3_bind_null (stmt_index, 1); assert (rc == SQLITE_OK);
            rc = sqlite3_bind_double (stmt_index, 2, no_value_present.f64); assert (rc == SQLITE_OK);
            break;

        case STITCH_NO_TYPE:
        default:
            fprintf (stderr, "Line: %d invalid field type found: %d\n", __LINE__, type);
            assert (0);
            break;
    }
    rc = sqlite3_bind_int64 (stmt_index, 3, field_id); assert (rc == SQLITE_OK);

    do
    {
        rc = sqlite3_step (stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY && rc != SQLITE_ROW && rc != SQLITE_DONE)
        {
            fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    rc = sqlite3_finalize (stmt_index);

cleanup:
    return rc;
}

int stitch_query_field (const StitchFile * file, const char * label, int64_t * field_id)
{
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;

    if (file->rank == 0)
    {
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select field_id from fields where field = ?", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
        rc = sqlite3_bind_text (stmt_index, 1, strdup (label), strlen (label), free); assert (rc == SQLITE_OK);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        if (rc == SQLITE_ROW)
        {
            //uint32_t db_type = 0;
            //int32_t db_length = 0;
            //const unsigned char * name = NULL;

            *field_id = sqlite3_column_int64 (stmt_index, 0);
            rc = sqlite3_finalize (stmt_index);
        }
        else
        {
            if (rc == SQLITE_DONE)
            {
                //printf ("found nothing\n");
                // close out the query before the insert
                rc = sqlite3_finalize (stmt_index);
                *field_id = STITCH_FIELD_NOT_FOUND_ID;
            }
            else
            {
                *field_id = STITCH_FIELD_NOT_FOUND_ID;
                // what's going on here?
                fprintf (stderr, "Line: %d rc: %d\n", __LINE__, rc);
                rc = sqlite3_finalize (stmt_index);
                sqlite3_close (file->db);
                assert (0);
            }
        }
        //printf ("rank: %d new field_id: %" PRId64 "\n", file->rank, *field_id);
    }
#ifdef STITCH_PARALLEL
    MPI_Bcast (field_id, 1, MPI_INT64_T, 0, file->comm);
#endif

cleanup:

    return rc;
}

static int stitch_write_block (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const void * state, int32_t * is_new_time_flag);
int stitch_write_block_int32 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const int32_t * state, int32_t * is_new_time_flag)
{
    return stitch_write_block (file, field_id, time, bb, (void *) state, is_new_time_flag);
}

int stitch_write_block_int64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const int64_t * state, int32_t * is_new_time_flag)
{
    return stitch_write_block (file, field_id, time, bb, (void *) state, is_new_time_flag);
}

int stitch_write_block_float64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const double * state, int32_t * is_new_time_flag)
{
    return stitch_write_block (file, field_id, time, bb, (void *) state, is_new_time_flag);
}

static int stitch_write_block (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const void * state, int32_t * is_new_time_flag)
{
    assert (file);
    assert (time);
    assert (bb);
    assert (state);
    assert (is_new_time_flag);
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
//    char * ErrMsg = NULL;
    char * state_buffer = NULL;
    int32_t x1 = bb [0];
    int32_t x2 = bb [1];
    int32_t y1 = bb [2];
    int32_t y2 = bb [3];
    int32_t z1 = bb [4];
    int32_t z2 = bb [5];
    //printf ("x1: %d x2: %d y1: %d y2: %d z1: %d z2: %d\n", x1, x2, y1, y2, z1, z2);

    int64_t timestamp = 0;
    double real_time = *time;
    int32_t new_time = false;

    if (field_id == STITCH_FIELD_NOT_FOUND_ID)
    {
        fprintf (stderr, "stitch_write_block: Field ID is invalid: %" PRId64 "\n", field_id);
        goto cleanup;
    }

    //printf ("%-6.1lf: (%lld, %lld, %lld) - (%lld, %lld, %lld)\n", timestep, x1, y1, z1, x2, y2, z2);
    //printf ("file: %lld state [0]: %lld [1]: %lld\n", (int) file, state [0], state [1]);
    rc = stitch_create_timestamp (file, &real_time, &timestamp, &new_time);
    *is_new_time_flag = new_time;
    //printf ("write_block: %3.1g, %lld, %d\n", real_time, timestamp, new_time == true);

    // we don't want the last value for each one which is why we don't do + 1
    int32_t state_size = ((x2 - x1) * (y2 - y1) * (z2 - z1));
#if defined(STITCH_PARALLEL)
    MPI_Status status;
    if (file->rank != 0)
    {
        MPI_Recv (&state_size, 1, MPI_INT, file->rank - 1, STITCH_MPI_TAG_WRITE_BLOCK, file->comm, &status);
    }
    else // get the size of the field and that is what we send as our value in the serialization
    {
#endif
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select type, length from fields where field_id = ?", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        rc = sqlite3_bind_int64 (stmt_index, 1, field_id);

        do
        {
            rc = sqlite3_step (stmt_index);
            //printf ("in da loop rc: %d\n", rc);
            if (rc != SQLITE_OK && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

	if (rc == SQLITE_DONE)
        {
            fprintf (stderr, "Field ID not found: %" PRId64 ". Line: %d rc = %d\n", field_id, __LINE__, rc);
            sqlite3_finalize (stmt_index);
            goto cleanup;
        }

        uint32_t field_type = STITCH_NO_TYPE;
        uint32_t field_length = 1;

        field_type = sqlite3_column_int (stmt_index, 0);
        field_length = sqlite3_column_int (stmt_index, 1);

        sqlite3_finalize (stmt_index);

        switch (field_type)
        {
            case STITCH_INT32:
                state_size *= 4;
                break;

            case STITCH_INT64:
            case STITCH_FLOAT64:
                state_size *= 8;
                break;

            case STITCH_NO_TYPE:
                fprintf (stderr, "Line: %d field_id: %" PRId64 " invalid field type found: %d\n", __LINE__, field_id, field_type);
                break;
        }
        assert (field_length >= 0);
        state_size *= field_length; // this accounts for all of the bytes in the type and the number of entries in the tensor
        //printf ("write field_id: %" PRId64 " type: %d field_length: %d state_size: %d\n", field_id, field_type, field_length, state_size);
#if defined(STITCH_PARALLEL)
    }
#endif

    do
    {
        rc = sqlite3_prepare_v2 (file->db, "insert into blocks (timestamp, field_id, x_min, y_min, z_min, x_max, y_max, z_max, state) values (?, ?, ?, ?, ?, ?, ?, ?, ?)", -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

    state_buffer = malloc (state_size);
    memcpy (state_buffer, state, state_size);

    rc = sqlite3_bind_int (stmt_index, 1, timestamp); if (rc != SQLITE_OK) fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db)); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int64 (stmt_index, 2, field_id); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 3, x1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 4, y1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 5, z1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 6, x2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 7, y2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 8, z2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_blob (stmt_index, 9, state_buffer, state_size, free); assert (rc == SQLITE_OK);

    do
    {
        rc = sqlite3_step (stmt_index);
        //printf ("in da loop rc: %d\n", rc);
        if (rc != SQLITE_OK && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

    rc = sqlite3_finalize (stmt_index);
//    int rowid = sqlite3_last_insert_rowid (file->db);
//printf ("rowid for chunk: %d\n", rowid);

#if defined(STITCH_PARALLEL)
    if (file->rank < file->size - 1)
    {
        MPI_Send (&state_size, 1, MPI_INT, file->rank + 1, STITCH_MPI_TAG_WRITE_BLOCK, file->comm);
    }
#endif

#ifdef STITCH_PARALLEL
    MPI_Barrier (file->comm);  // since we may read part of the group write, need to wait for everyone
#endif

cleanup:

    return rc;
}

#if 0
static void print_block (const char * prefix, int32_t * bb, int32_t * state)
{
    int32_t x1 = bb [0];
    int32_t y1 = bb [1];
    int32_t z1 = bb [2];
    int32_t x2 = bb [3];
    int32_t y2 = bb [4];
    int32_t z2 = bb [5];
    int32_t offset = 0;
    printf ("%s ", prefix);
    printf ("(%d, %d, %d) - (%d, %d, %d)\n", x1, y1, z1, x2, y2, z2);
    for (int32_t z = z1; z < z2; z++)
    {
        printf ("z(%d) ", z);
        for (int32_t y = y1; y < y2; y++)
        {
            printf ("y(%d) ", y);
            for (int32_t x = x1; x < x2; x++)
            {
                printf ("%d ", state [offset++]);
            }
            printf ("\n");
        }
    }
}
#endif

static int stitch_read_block (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, void * buffer, int32_t * is_new_time_flag);
int stitch_read_block_int32 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, int32_t * buffer, int32_t * is_new_time_flag)
{
    return stitch_read_block (file, field_id, time, bb, (void *) buffer, is_new_time_flag);
}

int stitch_read_block_int64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, int64_t * buffer, int32_t * is_new_time_flag)
{
    return stitch_read_block (file, field_id, time, bb, (void *) buffer, is_new_time_flag);
}

int stitch_read_block_float64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, double * buffer, int32_t * is_new_time_flag)
{
    return stitch_read_block (file, field_id, time, bb, (void *) buffer, is_new_time_flag);
}

static int stitch_read_block (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, void * buffer, int32_t * is_new_time_flag)
{
    assert (file);
    assert (time);
    assert (bb);
    assert (buffer);
    assert (is_new_time_flag);
    int rc = 0;
    sqlite3_stmt * stmt_index = NULL;
    const char * tail_index = NULL;
//    char * ErrMsg = NULL;
    int32_t x1 = bb [0];
    int32_t x2 = bb [1];
    int32_t y1 = bb [2];
    int32_t y2 = bb [3];
    int32_t z1 = bb [4];
    int32_t z2 = bb [5];
    void * new_block = 0;

    int64_t timestamp = -1;
    double real_time = *time;
    int32_t new_time = false;
    *is_new_time_flag = new_time;

    // we don't want the last value for each one which is why we don't do + 1
    int32_t size_to_nuke = (x2 - x1) * (y2 - y1) * (z2 - z1);
    int32_t field_info [2] = {STITCH_NO_TYPE, 1};
    int32_t no_value_i32 = STITCH_NO_VALUE;
    int64_t no_value_i64 = STITCH_NO_VALUE;
    double no_value_f64 = STITCH_NO_VALUE;
#ifdef STITCH_TIMING
    int max_timing = 2000;
    struct timeval time_start, time_end, time_diff;
    double timing [max_timing];
    int timing_offset = 0;
#endif
    // these 3 vars are used to get which times are wanted to avoid a join
    int32_t number_of_times = 0;
    int32_t times_selected = 0;
    int32_t * times = NULL;
#ifdef STITCH_PARALLEL
    if (file->rank == 0)
    {
#endif
#ifdef STITCH_TIMING
        gettimeofday (&time_start, NULL);
#endif
        do
        {
            rc = sqlite3_prepare_v2 (file->db, "select type, length, no_value_int, no_value_real from fields where field_id = ?", -1, &stmt_index, &tail_index);
            if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
#ifdef STITCH_TIMING
        gettimeofday (&time_end, NULL);
        timing [timing_offset++] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;

        gettimeofday (&time_start, NULL);
#endif
        rc = sqlite3_bind_int64 (stmt_index, 1, field_id);

        do
        {
            rc = sqlite3_step (stmt_index);
            //printf ("in da loop rc: %d\n", rc);
            if (rc != SQLITE_OK && rc != SQLITE_DONE && rc != SQLITE_ROW && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

        field_info [0] = sqlite3_column_int (stmt_index, 0);
        field_info [1] = sqlite3_column_int (stmt_index, 1);
        switch (field_info [0])
        {
            case STITCH_INT32: // we are relying on little-endian-ness. Pretty much everything does that today
                no_value_i32 = sqlite3_column_int (stmt_index, 2);
                break;

            case STITCH_INT64:
                no_value_i64 = sqlite3_column_int64 (stmt_index, 2);
                break;

            case STITCH_FLOAT64:
                no_value_f64 = sqlite3_column_double (stmt_index, 3);
                break;

            case STITCH_NO_VALUE:
            default:
                fprintf (stderr, "Line: %d invalid field type found: %d\n", __LINE__, field_info [0]);
                assert (0);
                break;
        }

        sqlite3_finalize (stmt_index);
#ifdef STITCH_TIMING
        gettimeofday (&time_end, NULL);
        timing [timing_offset++] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;
#endif
#ifdef STITCH_PARALLEL
    }
    MPI_Bcast (field_info, 2, MPI_INT, 0, file->comm);
    MPI_Bcast (&no_value_i32, 1, MPI_INT, 0, file->comm);
    MPI_Bcast (&no_value_i64, 1, MPI_INT64_T, 0, file->comm);
    MPI_Bcast (&no_value_f64, 1, MPI_DOUBLE, 0, file->comm);
    // do not need to serialize since reading does not cause locks
#endif

    int32_t * buffer_i32 = (int32_t *) buffer;
    int64_t * buffer_i64 = (int64_t *) buffer;
    double * buffer_f64 = (double *) buffer;

    //printf ("size to nuke: %d\n", size_to_nuke);
    switch (field_info [0])
    {
        case STITCH_INT32:
            for (int32_t i = 0; i < size_to_nuke; i++)
                buffer_i32 [i] = no_value_i32;
            break;

        case STITCH_INT64:
            for (int32_t i = 0; i < size_to_nuke; i++)
                buffer_i64 [i] = no_value_i64;
            break;

        case STITCH_FLOAT64:
            for (int32_t i = 0; i < size_to_nuke; i++)
                buffer_f64 [i] = no_value_f64;
            break;


        default:
            fprintf (stderr, "Line: %d invalid field type found: %d\n", __LINE__, field_info [0]);
            assert (0);
            break;
    }

#ifdef STITCH_TIMING
    gettimeofday (&time_start, NULL);
#endif
    rc = stitch_get_timestamp (file, &real_time, &timestamp, &new_time);
#ifdef STITCH_TIMING
    gettimeofday (&time_end, NULL);
    timing [timing_offset++] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;
#endif
    *is_new_time_flag = new_time;

// do I get the list in descending or ascending order?
// if descending, just copy the relevant pieces into the output buffer, but
// more copies. If ascending, then copy the whole last piece and then slice
// in just the pieces that are needed from each new block careful to not
// overwrite the existing data. Descending seems easier for now. Performance
// requirements may change this. I don't think I need the timestep if I do it
// in ascending order. Need to test this.
#ifdef STITCH_TIMING
    gettimeofday (&time_start, NULL);
#endif
    do
    {
        const char * query = NULL;
        query = "select count (*) from times";
        rc = sqlite3_prepare_v2 (file->db, query, -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
    do
    {
        rc = sqlite3_step (stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_ROW && rc != SQLITE_DONE)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    number_of_times = sqlite3_column_int (stmt_index, 0);
    rc = sqlite3_finalize (stmt_index);
    times = malloc (sizeof (int32_t) * number_of_times);

    do
    {
        const char * query = NULL;
        query = "select timestamp from times, globals where times.time < globals.absolute_tolerance + ?";
        rc = sqlite3_prepare_v2 (file->db, query, -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
    rc = sqlite3_bind_double (stmt_index, 1, *time); assert (rc == SQLITE_OK);
    do
    {
        rc = sqlite3_step (stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_ROW && rc != SQLITE_DONE)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    while (rc == SQLITE_ROW)
    {
        times [times_selected++] = sqlite3_column_int (stmt_index, 0);
        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    }
    rc = sqlite3_finalize (stmt_index);
    // now that we have a list of matching timestamps, we can add them to avoid the join.
    // Optimization: If the number of items == the total number, we can just eliminate that clause.

    do
    {
        // AABB (Axis Aligned Bounding Box) Intersection
        // you must check all axes INDEPENDENTLY.
        // https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
#ifndef SQLITE_MAX_SQL_LENGTH
#define SQLITE_MAX_SQL_LENGTH 1000000
#endif
#define SQLITE_MAX_SQL_LENGTH_CAT 998000
        char query [SQLITE_MAX_SQL_LENGTH] = "";
        const char * q_txt = 
                "select timestamp, field_id, x_min, y_min, z_min, x_max, y_max, z_max, state "
                "from blocks "
                "where "
                "    x_min <= ? and x_max >= ? "
                "and y_min <= ? and y_max >= ? "
                "and z_min <= ? and z_max >= ? "
                "and field_id = ? ";
        strncpy (query, q_txt, SQLITE_MAX_SQL_LENGTH_CAT);
        if (number_of_times != times_selected)
        {
            strncat (query, "and timestamp in (?", SQLITE_MAX_SQL_LENGTH_CAT);
            for (int i = 1; i < times_selected; i++)
            {
                strncat (query, ",?", SQLITE_MAX_SQL_LENGTH_CAT);
            }
            strncat (query, ") ", SQLITE_MAX_SQL_LENGTH_CAT);
        }
        else
        {
            times_selected = 0; // set to skip binding these parameters
        }
        strncat (query, "order by timestamp ", SQLITE_MAX_SQL_LENGTH_CAT);

        rc = sqlite3_prepare_v2 (file->db, query, -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_BUSY || rc == SQLITE_LOCKED);
#ifdef STITCH_TIMING
    gettimeofday (&time_end, NULL);
    timing [timing_offset++] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;
#endif

    rc = sqlite3_bind_int (stmt_index, 1, x2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 2, x1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 3, y2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 4, y1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 5, z2); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int (stmt_index, 6, z1); assert (rc == SQLITE_OK);
    rc = sqlite3_bind_int64 (stmt_index, 7, field_id); assert (rc == SQLITE_OK);
    for (int i = 0; i < times_selected; i++)
    rc = sqlite3_bind_int (stmt_index, 8 + i, times [i]); assert (rc == SQLITE_OK);
    //rc = sqlite3_bind_double (stmt_index, 8, *time); assert (rc == SQLITE_OK);
//printf ("time: %g bounding box x1: %d y1: %d z1: %d x2: %d y2: %d z2: %d\n", *time, x1, y1, z1, x2, y2, z2);
//if (file->rank == 0) printf ("1: %d 2: %d 3: %d 4: %d 5: %d 6: %d 7: %lld 8: %f\n", x2, x1, y2, y1, z2, z1, field_id, *time);
#ifdef STITCH_TIMING
    gettimeofday (&time_start, NULL);
#endif
    do
    {
        rc = sqlite3_step (stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_ROW && rc != SQLITE_DONE)
        {
            fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
            sqlite3_close (file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
#ifdef STITCH_TIMING
    gettimeofday (&time_end, NULL);
    timing [timing_offset++] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;
#endif

    //printf ("rank: %d timestamp: %" PIRd64 " rc: %d\n", file->rank, timestamp, rc);

    int32_t * new_block_i32 = 0;
    int64_t * new_block_i64 = 0;
    double * new_block_f64 = 0;
//int32_t min_val = 2000000000;
//int32_t max_val = -1;
#ifdef STITCH_TIMING
    gettimeofday (&time_start, NULL);
    int time_saved = timing_offset++;
#endif
    while (rc == SQLITE_ROW)
    {
        //double block_time = 0.0;
        int32_t new_block_x1 = 0;
        int32_t new_block_y1 = 0;
        int32_t new_block_z1 = 0;
        int32_t new_block_x2 = 0;
        int32_t new_block_y2 = 0;
        int32_t new_block_z2 = 0;

#ifdef STITCH_TIMING
        struct timeval start1, end1;
        gettimeofday (&time_start, NULL);
#endif
        //block_time = sqlite3_column_double (stmt_index, 0);
        field_id = sqlite3_column_int64 (stmt_index, 1);
        new_block_x1 = sqlite3_column_int (stmt_index, 2);
        new_block_y1 = sqlite3_column_int (stmt_index, 3);
        new_block_z1 = sqlite3_column_int (stmt_index, 4);
        new_block_x2 = sqlite3_column_int (stmt_index, 5);
        new_block_y2 = sqlite3_column_int (stmt_index, 6);
        new_block_z2 = sqlite3_column_int (stmt_index, 7);
#ifdef STITCH_TIMING
        gettimeofday (&time_end, NULL);
        timing [timing_offset++] = (end1.tv_sec - start1.tv_sec) + ((double) end1.tv_usec - start1.tv_usec)/1000000;
#endif

//if (file->rank == 0) printf ("new block (x1: %d, y1: %d, z1: %d) - (x2: %d, y2: %d, z2: %d)\n", new_block_x1, new_block_y1, new_block_z1, new_block_x2, new_block_y2, new_block_z2);

//        printf ("block_time: %g time: %g real_time: %g\n", block_time, *time, real_time);

        // proper way to read a blob
        uint32_t block_size = 0;
#ifdef STITCH_TIMING
        gettimeofday (&start1, NULL);
#endif
        block_size = sqlite3_column_bytes (stmt_index, 8);
        //printf ("block_size: %u\n", block_size);
        new_block = realloc (new_block, block_size); assert (new_block);
#ifdef STITCH_TIMING
        gettimeofday (&end1, NULL);
        timing [timing_offset++] = (end1.tv_sec - start1.tv_sec) + ((double) end1.tv_usec - start1.tv_usec)/1000000;
#endif
        if (!new_block)
        {
            goto cleanup;
        }
        memcpy (new_block, sqlite3_column_blob (stmt_index, 8), block_size);
        new_block_i32 = (int32_t *) new_block;
        new_block_i64 = (int64_t *) new_block;
        new_block_f64 = (double *) new_block;

        // do the copy of the right pieces into the 'buffer' parameter
        int32_t new_block_x_pos = x1;
        int32_t new_block_y_pos = y1;
        int32_t new_block_z_pos = z1;
        uint32_t new_block_offset = 0;
        uint32_t buffer_offset = 0;
#if 0
        {
        int bb1 [6];

        bb1 [0] = new_block_x1;
        bb1 [1] = new_block_y1;
        bb1 [2] = new_block_z1;
        bb1 [3] = new_block_x2;
        bb1 [4] = new_block_y2;
        bb1 [5] = new_block_z2;
        if (x1 == 0 && y1 == 86 && z1 == 0 && x2 == 127 && y2 == 344 && z2 == 30)
        {
            //print_block ("read_block", bb1, new_block);
            printf ("x1: %d y1: %d z1: %d new_block_x2: %d new_block_y2: %d new_block_z2: %d\n", x1, y1, z1, new_block_x1, new_block_y1, new_block_z1);
        }
        }
#endif

#define POINT_IN_BUFFER(x1,y1,z1,x2,y2,z2,x3,y3,z3) (x3 >= x1 ? x3 < x2 ? y3 >= y1 ? y3 < y2 ? z3 >= z1 ? z3 < z2 ? true : false : false : false : false : false: false)
#ifdef STITCH_TIMING
        gettimeofday (&start1, NULL);
#endif
        for (new_block_z_pos = new_block_z1; new_block_z_pos < new_block_z2; new_block_z_pos++)
        {
            for (new_block_y_pos = new_block_y1; new_block_y_pos < new_block_y2; new_block_y_pos++)
            {
                for (new_block_x_pos = new_block_x1; new_block_x_pos < new_block_x2; new_block_x_pos++)
                {
#if 0
                    if (x1 == 0 && y1 == 86 && z1 == 0 && x2 == 127 && y2 == 344 && z2 == 30)
printf ("POINT_IN_BUFFER (%d,%d,%d)-(%d,%d,%d) (%d,%d,%d) ", x1,y1,z1,x2,y2,z2,new_block_x_pos,new_block_y_pos, new_block_z_pos);
#endif
                    if (POINT_IN_BUFFER (x1,y1,z1,x2,y2,z2,new_block_x_pos,new_block_y_pos, new_block_z_pos))
                    {
#if 0
printf ("y\n");
#endif
                        // we don't want the last value for each one
                        // which is why we don't do + 1
                        buffer_offset = (new_block_x_pos - x1)
                                      + (new_block_y_pos - y1) * (x2 - x1)
                                      + (new_block_z_pos - z1) * (y2 - y1) * (x2 - x1);
                        buffer_offset *= field_info [1];
#if 0
                    if (x1 == 0 && y1 == 86 && z1 == 0 && x2 == 127 && y2 == 344 && z2 == 30)
printf ("buffer_offset = %d new_block_offset: %d read_val: %d\n", buffer_offset, new_block_offset, new_block [new_block_offset]);
#endif
                        assert (buffer_offset < size_to_nuke);
                        assert (new_block_offset < block_size);
                        switch (field_info [0])
                        {
                            case STITCH_INT32:
                                for (int32_t i = 0; i < field_info [1]; i++)
                                    buffer_i32 [buffer_offset + i] = new_block_i32 [new_block_offset + i];
                                break;

                            case STITCH_INT64:
                                for (int32_t i = 0; i < field_info [1]; i++)
                                    buffer_i64 [buffer_offset + i] = new_block_i64 [new_block_offset + i];
                                break;

                            case STITCH_FLOAT64:
                                for (int32_t i = 0; i < field_info [1]; i++)
                                    buffer_f64 [buffer_offset + i] = new_block_f64 [new_block_offset + i];
                                break;
                        }
//min_val = (new_block [new_block_offset] < min_val ? new_block [new_block_offset] : min_val);
//max_val = (new_block [new_block_offset] > max_val ? new_block [new_block_offset] : max_val);
                    }
#if 0
else
{
printf ("n\n");
}
#endif
                    new_block_offset += field_info [1];
                }
            }
        }
#ifdef STITCH_TIMING
        gettimeofday (&end1, NULL);
        timing [timing_offset++] = (end1.tv_sec - start1.tv_sec) + ((double) end1.tv_usec - start1.tv_usec)/1000000;
#endif

        do
        {
            rc = sqlite3_step (stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_ROW && rc != SQLITE_DONE)
            {
                fprintf (stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg (file->db));
                sqlite3_close (file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    }
#ifdef STITCH_TIMING
    gettimeofday (&time_end, NULL);
    timing [time_saved] = (time_end.tv_sec - time_start.tv_sec) + ((double) time_end.tv_usec - time_start.tv_usec)/1000000;
#endif

#if 0
printf ("copied min: %d max: %d\n", min_val, max_val);
min_val = 2000000000;
max_val = -100;
for (int i = 0; i < (x2 - x1) * (y2 - y1) * (z2 - z1); i++)
{
min_val = (buffer [i] < min_val ? buffer [i] : min_val);
max_val = (buffer [i] > max_val ? buffer [i] : max_val);
}
printf ("in array min: %d max: %d\n", min_val, max_val);
#endif

    rc = sqlite3_finalize (stmt_index);

cleanup:
    free (new_block);
    free (times);
#ifdef STITCH_TIMING
    if (file->rank == 0)
    {
        for (int i = 0; i < timing_offset; i++)
            printf ("time %d: %lf\n", i, timing [i]);
    }
#endif

    return 0;
}

int stitch_get_global_bounds(const StitchFile * file, int64_t field_id, int32_t * bb)
{
    assert(file);
    assert(bb);
    /*
        bb [0] = x1;
        bb [1] = y1;
        bb [2] = z1;
        bb [3] = x2;
        bb [4] = y2;
        bb [5] = z2;

    */
    int rc = 0;
    sqlite3_stmt* stmt_index = NULL;
    const char* tail_index = NULL;

    do
    {
        rc = sqlite3_prepare_v2(file->db, "select min(blocks.x_min), max(blocks.x_max), min(blocks.y_min), max(blocks.y_max), min(blocks.z_min), max(blocks.z_max) from blocks where field_id = ?", -1, &stmt_index, &tail_index);
        if (rc != SQLITE_OK && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf(stderr, "%d get_global_bounds, prepare Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg(file->db));
            sqlite3_close(file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    rc = sqlite3_bind_int64 (stmt_index, 1, field_id); assert (rc == SQLITE_OK);
    do
    {
        rc = sqlite3_step(stmt_index);
        if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE && rc != SQLITE_LOCKED && rc != SQLITE_BUSY)
        {
            fprintf(stderr, "%d get_global_bounds, prepare Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg(file->db));
            sqlite3_close(file->db);
            goto cleanup;
        }
    } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);

    while (rc == SQLITE_ROW)
    {
        bb [0] = sqlite3_column_int(stmt_index, 0);
        bb [1] = sqlite3_column_int(stmt_index, 1);
        bb [2] = sqlite3_column_int(stmt_index, 2);
        bb [3] = sqlite3_column_int(stmt_index, 3);
        bb [4] = sqlite3_column_int(stmt_index, 4);
        bb [5] = sqlite3_column_int(stmt_index, 5);
        do
        {
            rc = sqlite3_step(stmt_index);
            if (rc != SQLITE_OK && rc != SQLITE_BUSY && rc != SQLITE_LOCKED && rc != SQLITE_DONE && rc != SQLITE_ROW)
            {
                fprintf(stderr, "Line: %d SQL error (%d): %s\n", __LINE__, rc, sqlite3_errmsg(file->db));
                sqlite3_close(file->db);
                goto cleanup;
            }
        } while (rc == SQLITE_LOCKED || rc == SQLITE_BUSY);
    }

cleanup:
    rc = sqlite3_finalize(stmt_index);
    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "%d Line: %d SQL error (%d): %s\n", file->rank, __LINE__, rc, sqlite3_errmsg(file->db));
        sqlite3_close(file->db);
    }

    return 0;
}
