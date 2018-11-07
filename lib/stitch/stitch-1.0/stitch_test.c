// moving parallel will be step 2
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "stitch.h"

void print_block (const char * prefix, int32_t * bb, int32_t * state)
{
    int32_t x1 = bb [0];
    int32_t x2 = bb [1];
    int32_t y1 = bb [2];
    int32_t y2 = bb [3];
    int32_t z1 = bb [4];
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

int main (int argc, char ** argv)
{
    StitchFile * file;

    double time = 1.0; // goes up by 1 each time
    double time_increment = 1.0;

    int increment = 2; // how much to move for each block

    int32_t g_x1 = 1; // global size
    int32_t g_y1 = 1;
    int32_t g_z1 = 1;
    int32_t g_x2 = 8; // global size
    int32_t g_y2 = 8;
    int32_t g_z2 = 8;

    int block_size = 4; // how big to get in each dimension

    int32_t x1 = 0;
    int32_t y1 = 0;
    int32_t z1 = 0;
    int32_t x2 = 0;
    int32_t y2 = 0;
    int32_t z2 = 0;

    const char * field = "spin";
    int64_t field_id = 0;

    // 3-D = 3 x block_size
    int32_t * state = (int32_t *) malloc (sizeof (int32_t) * block_size * block_size * block_size);

    MPI_Init (&argc, &argv);

    stitch_open (&file, MPI_COMM_WORLD, "test.st");
//    stitch_open (&file, "test.st");

    stitch_set_parameters (file, 1.0e-9, 1.0e-15, -1);
    union StitchTypesUnion no_value;
    no_value.i32 = STITCH_NO_VALUE;

    stitch_create_field (file, field, STITCH_INT32, no_value, 1, &field_id);
    assert (field_id != STITCH_FIELD_NOT_FOUND_ID);

    // just to test, query it again
    stitch_query_field (file, field, &field_id);
    assert (field_id != STITCH_FIELD_NOT_FOUND_ID);

    for (int32_t z = g_z1; z < g_z2; z += increment)
    {
        for (int32_t y = g_y1; y < g_y2; y += increment)
        {
            for (int32_t x = g_x1; x < g_x2; x += increment)
            {
                double absolute_tolerance;
                double relative_tolerance;
                int64_t no_value_present;
                double first_time;
                double last_time;
                int32_t bb [6];
                int32_t new_time = 0;

                stitch_get_parameters (file, &absolute_tolerance, &relative_tolerance, &no_value_present, &first_time, &last_time);
//                printf ("last_time: %-6.1lf\n", last_time);

                bb [0] = x;
                bb [1] = x + block_size;
                bb [2] = y;
                bb [3] = y + block_size;
                bb [4] = z;
                bb [5] = z + block_size;
                if (bb [1] > g_x2)
                    bb [1] = g_x2;
                if (bb [3] > g_y2)
                    bb [3] = g_y2;
                if (bb [5] > g_z2)
                    bb [5] = g_z2;
printf ("read: (%d, %d, %d) - (%d, %d, %d)\n", bb [0], bb [2], bb [4], bb [1], bb [3], bb [5]);
                stitch_read_block_int32 (file, field_id, &time, bb, state, &new_time);

//                print_block ("read", bb, state);

                // check the state for -1 values saying that it needs init

                // do compute
                int32_t offset = 0;
                int32_t new_val = (int32_t) time;
//printf ("new_val: %d\n", new_val);
                for (int32_t n_z = z; n_z < z + block_size; n_z++)
                {
                    for (int32_t n_y = y; n_y < y + block_size; n_y++)
                    {
                        for (int32_t n_x = x; n_x < x + block_size; n_x++)
                        {
//printf ("offset: %d\n", offset);
                            state [offset++] = new_val;
                        }
                    }
                }
//printf ("offset: %d\n", offset);

               // write new block
               stitch_write_block_int32 (file, field_id, &time, bb, state, &new_time);
               time += time_increment;
            }
        }
    }

    stitch_close (&file);

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize ();

    return 0;
}
