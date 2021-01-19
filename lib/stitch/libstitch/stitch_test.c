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
// moving parallel will be step 2
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>

#include "stitch.h"

// to generate testing data, set this to 1. Otherwise output goes to stdout
#define STITCH_TESTING_DATA 0

void print_block (FILE * f, const char * prefix, int32_t * bb, int32_t * state)
{
    int32_t x1 = bb [0];
    int32_t x2 = bb [1];
    int32_t y1 = bb [2];
    int32_t y2 = bb [3];
    int32_t z1 = bb [4];
    int32_t z2 = bb [5];
    int32_t offset = 0;
#if STITCH_TESTING_DATA
    // when generating new test data, use this version of the function to get
    // the values after a proper visual inspection
    for (int32_t z = z1; z < z2; z++)
    {
        for (int32_t y = y1; y < y2; y++)
        {
            for (int32_t x = x1; x < x2; x++)
            {
                fprintf (f, "%2d,", state [offset++]);
            }
        }
    }
    fprintf (f,"\n");
#else
    fprintf (f, "%s ", prefix);
    fprintf (f, "(%d, %d, %d) - (%d, %d, %d)\n", x1, y1, z1, x2, y2, z2); 
    for (int32_t z = z1; z < z2; z++)
    {
        fprintf (f, "z(%d)\n", z);
        for (int32_t y = y1; y < y2; y++)
        {
            fprintf (f, "  y(%d) ", y);
            for (int32_t x = x1; x < x2; x++)
            {
                fprintf (f, "%3d ", state [offset++]);
            }
            fprintf (f, "\n");
        }
    }
#endif
}

// compare the read block against the key
int test_stitch_block (int rank, int32_t * bb, int32_t * state, int32_t * test_data)
{
    int32_t x1 = bb [0];
    int32_t x2 = bb [1];
    int32_t y1 = bb [2];
    int32_t y2 = bb [3];
    int32_t z1 = bb [4];
    int32_t z2 = bb [5];
    int32_t offset = 0;
    for (int32_t z = z1; z < z2; z++)
    {
        for (int32_t y = y1; y < y2; y++)
        {
            for (int32_t x = x1; x < x2; x++)
            {
                if (state [offset] != test_data [offset])
                {
                    fprintf (stderr, "test failed. Rank: %d, (%d. %d, %d) %d does not match %d\n", rank, x, y, z, state [offset], test_data [offset]);
                    exit (-1);
                }
                offset++;
            }
        }
    }
    return 1;
}

int main (int argc, char ** argv)
{
    StitchFile * file;

    double time = 1.0; // goes up by 1 each time
    double time_increment = 1.0;

    int x_increment = 3; // how much to move for each block when we move in that direction
    int y_increment = 3; // how much to move for each block when we move in that direction
    int z_increment = 1; // how much to move for each block when we move in that direction

    int np_x = 2; // number of procs in x
    int np_y = 2; // number of procs in y
    int np_z = 1; // number of procs in z

    int32_t g_x1 = 0; // global size
    int32_t g_y1 = 0;
    int32_t g_z1 = 0;
    int32_t g_x2 = 13; // global size
    int32_t g_y2 = 13;
    int32_t g_z2 = 4; // do a single layer

    int x_block_size = 2; // how big to get in each dimension
    int y_block_size = 2; // how big to get in each dimension
    int z_block_size = 3; // how big to get in each dimension

    int32_t x1 = 0;
    int32_t y1 = 0;
    int32_t z1 = 0;
    int32_t x2 = 0;
    int32_t y2 = 0;
    int32_t z2 = 0;

#include "stitch_test_data.h"

    const char * field = "spin";
    int64_t field_id = 0;

    int32_t * state = (int32_t *) malloc (sizeof (int32_t) * x_block_size * y_block_size * z_block_size);
    int32_t * full_state = malloc ((g_x2 - g_x1) * (g_y2 - g_y1) * (g_z2 - g_z1) * sizeof (int32_t));

    int rank;
    int size;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    assert (size == 4);

    if (rank == 0)
    {
        FILE * f = NULL;

#if STITCH_TESTING_DATA
        f = fopen ("final", "a");
#else
        f = stdout;
#endif
        fprintf (f, "nprocs in x, y, z: %d, %d, %d\n", np_x, np_y, np_z);
        fprintf (f, "block size in x, y, z: %d, %d, %d\n", x_block_size, y_block_size, z_block_size);
        fprintf (f, "increment in x, y, z: %d, %d, %d\n", x_increment, y_increment, z_increment);
        fprintf (f, "total domain bounds (min) - (max): (%d, %d, %d) - (%d, %d, %d)\n", g_x1, g_y1, g_z1, g_x2, g_y2, g_z2);
        fprintf (f, "sweeping 0->n in X, then Y, then Z\n");
        unlink ("test.st");

#if STITCH_TESTING_DATA
        fclose (f);
#endif
    }

    MPI_Barrier (MPI_COMM_WORLD);

    // these are used to step through the testing data
    int xy_step = 0;
    int z_step = 0;

    stitch_open (&file, MPI_COMM_WORLD, "test.st");

    stitch_set_parameters (file, 1.0e-9, 1.0e-15, -1);
    union StitchTypesUnion no_value;
    no_value.i32 = STITCH_NO_VALUE;

    stitch_create_field (file, field, STITCH_INT32, no_value, 1, &field_id);
    assert (field_id != STITCH_FIELD_NOT_FOUND_ID);

    // just to test, query it again
    stitch_query_field (file, field, &field_id);
    assert (field_id != STITCH_FIELD_NOT_FOUND_ID);

    int32_t bb [6];
    int32_t new_time = 0;
    double last_time;
    for (int32_t z = g_z1; z + z_increment < g_z2; z += z_increment, z_step++)
    {
        for (int32_t y = g_y1; y + y_increment < g_y2; y += y_increment)
        {
            for (int32_t x = g_x1; x + x_increment < g_x2; x += x_increment, xy_step++)
            {
                double absolute_tolerance;
                double relative_tolerance;
                int64_t no_value_present;
                double first_time;
                new_time = 0;

                stitch_get_parameters (file, &absolute_tolerance, &relative_tolerance, &no_value_present, &first_time, &last_time);

                bb [0] = x + (rank % np_x) * np_x;
                bb [1] = bb [0] + x_block_size;
                bb [2] = y + ((rank / np_x) % np_y) * np_y;
                bb [3] = bb [2] + y_block_size;
                bb [4] = z + (rank / (np_x * np_y)) * np_z;
                bb [5] = bb [4] + z_block_size;
                stitch_read_block_int32 (file, field_id, &time, bb, state, &new_time);

#if 1
                //if (rank == 0)
                {
                    char rank_str [20]; 
#if STITCH_TESTING_DATA
                    sprintf (rank_str, "rank_%d", rank);
                    FILE * f = fopen (rank_str, "a");

                    print_block (f, rank_str, bb, state);

                    fclose (f);
#else
                    sprintf (rank_str, "rank_%d ", rank);
                    int32_t * test_data = NULL;
                    //print_block (stdout, rank_str, bb, state);
                    switch (rank)
                    {
                        case 0: test_data = test_rank_0 [xy_step]; break;
                        case 1: test_data = test_rank_1 [xy_step]; break;
                        case 2: test_data = test_rank_2 [xy_step]; break;
                        case 3: test_data = test_rank_3 [xy_step]; break;
                    }
                    test_stitch_block (rank, bb, state, test_data);
#endif
                }
#endif


                // do compute
                int32_t offset = 0;
                int32_t new_val = (int32_t) time;
                for (int32_t n_z = z; n_z < z + z_block_size; n_z++)
                {
                    for (int32_t n_y = y; n_y < y + y_block_size; n_y++)
                    {
                        for (int32_t n_x = x; n_x < x + x_block_size; n_x++)
                        {
                            state [offset++] = new_val;
                        }
                    }
                }

               // write new block
//printf ("[%d] write: (%d, %d, %d) - (%d, %d, %d)\n", rank, bb [0], bb [2], bb [4], bb [1], bb [3], bb [5]);
               stitch_write_block_int32 (file, field_id, &time, bb, state, &new_time);
               time += time_increment;
            }
        }

        bb [0] = g_x1;
        bb [1] = g_x2;
        bb [2] = g_y1;
        bb [3] = g_y2;
        bb [4] = g_z1;
        bb [5] = g_z2;
        double read_time = last_time + 1;

        // parallel call has to be done be all
        stitch_read_block_int32 (file, field_id, &read_time, bb, full_state, &new_time);

        if (rank == 0)
        {
#if STITCH_TESTING_DATA
            FILE * f = fopen ("final", "a");
            print_block (f, "final", bb, full_state);
            fclose (f);
#else
            //print_block (stdout, "final", bb, full_state);
            test_stitch_block (rank, bb, full_state, test_final [z_step]);
#endif
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (rank == 0)
        printf ("tests passed\n");

cleanup:
    stitch_close (&file);

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize ();

    return 0;
}
