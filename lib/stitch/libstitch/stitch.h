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
#include <stdint.h>

///////////////////////////////////////////////////////////////////////////////
// This is setup to build both serial and parallel libraries. To build
// serial libraries, just use stock. To build parallel, be sure to define
// STITCH_PARALLEL
// and it will shift to include the MPI pieces.
///////////////////////////////////////////////////////////////////////////////
#ifndef STITCH_H_
#define STITCH_H_

#ifdef __cplusplus
extern "C"{
#endif 

// bit 0 = time value (0 or 1)
// bits 1 and 2 = data found (00, 01, or 10) 11 is currently an error
enum STITCH_FLAGS
{
     STITCH_EXISTING_TIME =       0    // 0
    ,STITCH_UNKNOWN_TIME =        1    // 1
    ,STITCH_NO_DATA_FOUND =       0    // 00x
    ,STITCH_PARTIAL_DATA_FOUND =  2    // 01x
    ,STITCH_ALL_DATA_FOUND =      4    // 10x
};

#define STITCH_NO_VALUE -1
#define STITCH_DEFAULT_ABSOLUTE_TOLERANCE 1.0e-9
#define STITCH_DEFAULT_RELATIVE_TOLERANCE 1.0e-15
#define STITCH_FIELD_NOT_FOUND_ID -1

typedef struct StitchFile_struct StitchFile;

enum STITCH_TYPES
{
     STITCH_NO_TYPE = 0
    ,STITCH_INT32   = 1
    ,STITCH_INT64   = 2
    ,STITCH_FLOAT64 = 3
};

union StitchTypesUnion
{
    int32_t i32;
    int64_t i64;
    double f64;
};

#ifdef STITCH_PARALLEL
int stitch_open (StitchFile ** file, MPI_Comm comm, const char * filename);
#else
int stitch_open (StitchFile ** file, const char * filename);
#endif
int stitch_close (StitchFile ** file);

int stitch_set_parameters (const StitchFile * file, double absolute_tolerance, double relative_tolerance, int64_t default_no_value_present);

int stitch_get_parameters (const StitchFile * file, double * absolute_tolerance, double * relative_tolerance, int64_t * default_no_value_present, double * first_time, double * last_time);

int stitch_get_times_count (const StitchFile * file, int64_t * count);

int stitch_get_times (const StitchFile * file, int64_t * count, double * times);

// retrieve the number of fields
int stitch_get_fields_count (const StitchFile * file, int64_t * count);

// count is an in/out parameter. For IN, it is the size of the buffers. For OUT, it is the number of entries.
// the labels are allocated with malloc and must be freed by the caller.
int stitch_get_fields (const StitchFile * file, int64_t * count, int64_t * field_ids, char ** labels, enum STITCH_TYPES * types, int32_t * per_site_lengths, union StitchTypesUnion * no_value_present_values);

// the Python interface needs to determine which of the type-safe C functions to call based on the type of the field being written.
int stitch_get_field_type (const StitchFile * file, int64_t field_id, enum STITCH_TYPES * type);

// return value is a ID that corresponds to this field in the file (field_id).
// currently only a per_site_length of 1 is supported, but other lengths will be incorporated.
int stitch_create_field (const StitchFile * file, const char * label, enum STITCH_TYPES type, union StitchTypesUnion no_value_present, int32_t per_site_length, int64_t * field_id);
int stitch_set_field_no_value_present (const StitchFile * file, int64_t field_id, union StitchTypesUnion no_value_present);

// if label is not found, the field_id is set to STITCH_FIELD_NOT_FOUNT_ID
int stitch_query_field (const StitchFile * file, const char * label, int64_t * field_id);

//For write:
//write_flag == True when new time step created
//write_flag == False when new time step NOT created
int stitch_write_block_int32 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const int32_t * buffer, int32_t * write_flag);
int stitch_write_block_int64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const int64_t * buffer, int32_t * write_flag);
int stitch_write_block_float64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const double * buffer, int32_t * write_flag);

//For read:
//read_flag == True when time does not exist in the file, but data may be returned based on older times
//read_flag == False when exists in the database. The actual data returned will vary based on what is selected. It could be nothing.
int stitch_read_block_int32 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, int32_t * buffer, int32_t * is_new_time);
int stitch_read_block_int64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, int64_t * buffer, int32_t * is_new_time);
int stitch_read_block_float64 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, double * buffer, int32_t * is_new_time);

// discover the global bounds (in 3D) written so far
int stitch_get_global_bounds(const StitchFile * file, int64_t field_id, int32_t * bb);

#ifdef __cplusplus
}
#endif

#endif
