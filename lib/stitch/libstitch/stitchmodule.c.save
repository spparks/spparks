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
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>
#ifdef STITCH_PARALLEL
#include <mpi4py/mpi4py.h>
#include "mpi.h"
#endif
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include "stitch.h"


//int stitch_open (StitchFile ** file, const char * filename);
static PyObject * stitch_open_wrapper (PyObject * self, PyObject * args)
{
    int64_t * file = 0;
    char * filename = 0;
#ifdef STITCH_PARALLEL
    PyObject * comm_obj = 0;
    MPI_Comm * comm;
#endif

    PyObject * ret = NULL;

    // parse arguments
#ifdef STITCH_PARALLEL
    if (!PyArg_ParseTuple (args, "Os", &comm_obj, &filename))
#else
    if (!PyArg_ParseTuple (args, "s", &filename))
#endif
    {
        return NULL;
    }

#ifdef STITCH_PARALLEL
    comm = PyMPIComm_Get (comm_obj);
#endif

    // run the actual function
#ifdef STITCH_PARALLEL
    int rc = stitch_open ((StitchFile **) &file, *comm, filename);
#else
    int rc = stitch_open ((StitchFile **) &file, filename);
#endif

    ret = Py_BuildValue ("iL", rc, file);

    return ret;
}

//int stitch_close (StitchFile ** file);
static PyObject * stitch_close_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "L", &file))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_close ((StitchFile **) &file);

    ret = PyLong_FromLong (rc);

    return ret;
}

//int stitch_set_parameters (const StitchFile * file, double time_tolerance);
static PyObject * stitch_set_parameters_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    double absolute_tolerance = 0.0;
    double relative_tolerance = 0.0;
    int64_t no_value_present = STITCH_NO_VALUE;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LddL", &file, &absolute_tolerance, &relative_tolerance, &no_value_present))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_set_parameters ((StitchFile *) file, absolute_tolerance, relative_tolerance, no_value_present);

    //ret = PyLong_FromLong (rc);
    ret = Py_BuildValue ("i", rc);

    return ret;
}

//int stitch_get_parameters (const StitchFile * file, double * time_tolerance, double * first_time, double * last_time);
static PyObject * stitch_get_parameters_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    double absolute_tolerance = 0;
    double relative_tolerance = 0;
    int64_t no_value_present = 0;
    double first_time = 0;
    double last_time = 0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "L", &file))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_get_parameters ((StitchFile *) file, &absolute_tolerance, &relative_tolerance, &no_value_present, &first_time, &last_time);

    //ret = PyLong_FromLong (rc);
    ret = Py_BuildValue ("iddLdd", rc, absolute_tolerance, relative_tolerance, no_value_present, first_time, last_time);

    return ret;
}

#if 0
//int stitch_create_timestamp (const StitchFile * file, double time, int64_t * timestamp);
static PyObject * stitch_create_timestamp_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    int64_t timestamp = 0;
    double my_time = 0.0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LD", &file, &my_time))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_create_timestamp ((StitchFile *) file, my_time, &timestamp);

    ret = Py_BuildValue ("iL", rc, timestamp);

    return ret;
}
#endif

//int stitch_get_timestamps (const StitchFile * file, int64_t * count, int64_t * timestamps, double * times);
static PyObject * stitch_get_times_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    double * times = 0;
    int64_t count = 0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "L", &file))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_get_times_count ((StitchFile *) file, &count);

    // Build py array from scratch with fortran layout
    PyObject * times_ret = 0;

    int32_t nd = 1;
    npy_intp nyp_new_dims [1] = {count};
    int flags = NPY_ARRAY_FARRAY;

    PyArray_Descr * times_descr = PyArray_DescrFromType (NPY_DOUBLE);
    times_ret = PyArray_NewFromDescr (&PyArray_Type, times_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
    times = (double *) PyArray_DATA ((PyArrayObject *) times_ret);

    //printf ("count: %lld\n", count);

    rc = stitch_get_times ((StitchFile *) file, &count, times);

    ret = Py_BuildValue ("iN", rc, times_ret);

    return ret;
}

// get fields
static PyObject * stitch_get_fields_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    int64_t count = 0;
    int64_t * field_ids = 0;
    int32_t * types = 0;
    int32_t * lengths = 0;
    union StitchTypesUnion * no_value_presents = 0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "L", &file))
    {
        return NULL;
    }

    // run the actual function
    int rc = stitch_get_fields_count ((StitchFile *) file, &count);

    // Build py array from scratch with fortran layout
    PyObject * field_ids_ret = 0;
    PyObject * labels_ret = 0;
    PyObject * types_ret = 0;
    PyObject * lengths_ret = 0;

    int32_t nd = 1;
    npy_intp nyp_new_dims [1] = {count};
    int flags = NPY_ARRAY_FARRAY;

    PyArray_Descr * field_ids_descr = PyArray_DescrFromType (NPY_INT64);
    field_ids_ret = PyArray_NewFromDescr (&PyArray_Type, field_ids_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
    field_ids = (int64_t *) PyArray_DATA ((PyArrayObject *) field_ids_ret);

    labels_ret = PyList_New (count);

    PyArray_Descr * types_descr = PyArray_DescrFromType (NPY_INT32);
    types_ret = PyArray_NewFromDescr (&PyArray_Type, types_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
    types = (int32_t *) PyArray_DATA ((PyArrayObject *) types_ret);

    PyArray_Descr * lengths_descr = PyArray_DescrFromType (NPY_INT32);
    lengths_ret = PyArray_NewFromDescr (&PyArray_Type, lengths_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
    lengths = (int32_t *) PyArray_DATA ((PyArrayObject *) lengths_ret);

    no_value_presents = malloc (sizeof (union StitchTypesUnion) * count);
    char ** labels_temp = (char **) malloc (sizeof (char *) * count);

    rc = stitch_get_fields ((StitchFile *) file, &count, field_ids, labels_temp, (enum STITCH_TYPES *) types, lengths, no_value_presents);

    for (int i = 0; i < count; i++)
    {
        PyObject * str = PyUnicode_FromString (labels_temp [i]);
        free (labels_temp [i]);
        PyList_SetItem (labels_ret, i, str);
    }
    free (labels_temp);

    PyObject * no_value_presents_list = PyList_New (count);
    assert (no_value_presents_list);
    for (int32_t i = 0; i < count; i++)
    {
        PyObject * item = 0;
        switch (types [i])
        {
            case STITCH_INT32:
                item = PyLong_FromLong (no_value_presents [i].i32);
                break;
                
            case STITCH_INT64:
                item = PyLong_FromLongLong (no_value_presents [i].i64);
                break;
                
            case STITCH_FLOAT64:
                item = PyFloat_FromDouble (no_value_presents [i].f64);
                break;
        }
        PyList_SetItem (no_value_presents_list, i, item);
    }

    ret = Py_BuildValue ("iNNNNN", rc, field_ids_ret, labels_ret, types_ret, lengths_ret, no_value_presents_list);

    return ret;
}

//int stitch_create_field (const StitchFile * file, const char * label, enum STITCH_TYPES type, union StitchTypesUnion no_value_present, int32_t per_site_length, int64_t * field_id);
static PyObject * stitch_create_field_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    const char * label = 0;
    int32_t type = STITCH_NO_TYPE;
    PyObject * no_value_present_obj = 0;
    union StitchTypesUnion no_value_present;
    int32_t per_site_length = 1;
    int64_t field_id = STITCH_FIELD_NOT_FOUND_ID;
    int rc = 0;

    PyObject * ret = 0;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LsiiO", &file, &label, &type, &per_site_length, &no_value_present_obj))
    {
        return NULL;
    }

    //printf ("%ld\n", Py_TYPE (no_value_present_obj)->tp_itemsize);
    // while this is promising, it doesn't work. It gives 4 for any size int and 0 for floats.
    if (PyLong_Check (no_value_present_obj))
    {
        if (type == STITCH_INT32)
        {
            no_value_present.i32 = PyLong_AsLong (no_value_present_obj);
            //printf ("no_value_present.i32: %d\n", no_value_present.i32);
        }
        else
        {
            if (type == STITCH_INT64)
            {
                no_value_present.i64 = PyLong_AsLongLong (no_value_present_obj);
            //printf ("no_value_present.i64: %lld\n", no_value_present.i64);
            }
        }
    }
    else
    {
        if (PyFloat_Check (no_value_present_obj))
        {
            no_value_present.f64 = PyFloat_AS_DOUBLE (no_value_present_obj);
            //printf ("no_value_present.f64: %g\n", no_value_present.f64);
        }
        else
        {
            fprintf (stderr, "invalid object type pass to create_field for the no_value_present value\n");
            assert (0);
        }
    }

    rc = stitch_create_field ((StitchFile *) file, label, type, no_value_present, per_site_length, &field_id);

    ret = Py_BuildValue ("iL", rc, field_id);

    return ret;
}

//int stitch_set_field_no_value_present (const StitchFile * file, int64_t field_id, union StitchTypesUnion no_value_present);
static PyObject * stitch_set_field_no_value_present_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    int64_t field_id = 0;
    PyObject * no_value_present_obj = 0;
    union StitchTypesUnion no_value_present;
    int rc = 0;

    PyObject * ret = 0;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LLO", &file, &field_id, &no_value_present_obj))
    {
        return NULL;
    }

    switch (Py_SIZE (no_value_present_obj))
    {
        case 4:
            if (PyLong_Check (no_value_present_obj))
                no_value_present.i32 = PyLong_AsLong (no_value_present_obj);
            else
                assert (0);
            break;

        case 8:
            if (PyLong_Check (no_value_present_obj))
                no_value_present.i64 = PyLong_AsLongLong (no_value_present_obj);
            else
                if (PyFloat_Check (no_value_present_obj))
                    no_value_present.f64 = PyFloat_AS_DOUBLE (no_value_present_obj);
                else
                    assert (0);
            break;

        default:
            fprintf (stderr, "Line: %d size: %ld invalid\n", __LINE__, Py_SIZE (no_value_present_obj));
            assert (0);
            break;
    }

    rc = stitch_set_field_no_value_present ((StitchFile *) file, field_id, no_value_present);

    ret = Py_BuildValue ("i", rc);

    return ret;
}

// query field
static PyObject * stitch_query_field_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    const char * name = 0;
    int64_t field_id = 0;
    int rc = 0;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "Ls", &file, &name))
    {
        return NULL;
    }

    rc = stitch_query_field ((StitchFile *) file, name, &field_id);

    //ret = PyLong_FromLong (rc);
    ret = Py_BuildValue ("iL", rc, field_id);

    return ret;
}

//int stitch_write_block (const StitchFile * file, double timestep, int32_t * dims, const int32_t * buffer);
static PyObject * stitch_write_block_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    double time = 0;
    int32_t new_time = 0;
    /*
    int x1 = 0;
    int y1 = 0;
    int z1 = 0;
    int x2 = 0;
    int y2 = 0;
    int z2 = 0;
    */
    PyObject * dims_parse = NULL;
    PyObject * state_parse = NULL;

    int32_t * passed_dims = NULL;
    int64_t field_id = 0;
    int rc = 0;

    int32_t * state_i32 = NULL;
    int64_t * state_i64 = NULL;
    double * state_f64 = NULL;

    PyObject * ret = NULL;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LLdOO", &file, &field_id, &time, &dims_parse, &state_parse))
    {
        return NULL;
    }

    PyArrayObject * dim_obj = (PyArrayObject *) PyArray_FROM_OTF (dims_parse ,NPY_INT32, NPY_ARRAY_FARRAY);
    passed_dims = (int32_t *) PyArray_DATA (dim_obj);
    /*
    x1 = passed_dims [0];
    x2 = passed_dims [1];
    y1 = passed_dims [2];
    y2 = passed_dims [3];
    z1 = passed_dims [4];
    z2 = passed_dims [5];
    printf("WRITE A: xlo,xhi,ylo,yhi,zlo,zhi = %d,%d,%d,%d,%d,%d\n",x1,x2,y1,y2,z1,z2);
    */
    PyArrayObject * state_obj = 0;
    enum STITCH_TYPES type = STITCH_NO_TYPE;
    rc = stitch_get_field_type ((const StitchFile *) file, field_id, &type);

    // run the actual function
    switch (type)
    {
        case STITCH_INT32:
            state_obj = (PyArrayObject *) PyArray_FROM_OTF (state_parse, NPY_INT32, NPY_ARRAY_FARRAY);
            state_i32 = (int32_t *) PyArray_DATA (state_obj);
            rc = stitch_write_block_int32 ((const StitchFile *) file, field_id, &time, passed_dims, state_i32, &new_time);
            break;

        case STITCH_INT64:
            state_obj = (PyArrayObject *) PyArray_FROM_OTF (state_parse, NPY_INT64, NPY_ARRAY_FARRAY);
            state_i64 = (int64_t *) PyArray_DATA (state_obj);
            rc = stitch_write_block_int64 ((const StitchFile *) file, field_id, &time, passed_dims, state_i64, &new_time);
            break;

        case STITCH_FLOAT64:
            state_obj = (PyArrayObject *) PyArray_FROM_OTF (state_parse, NPY_FLOAT64, NPY_ARRAY_FARRAY);
            state_f64 = (double *) PyArray_DATA (state_obj);
            rc = stitch_write_block_float64 ((const StitchFile *) file, field_id, &time, passed_dims, state_f64, &new_time);
            break;

        case STITCH_NO_TYPE:
        default:
            fprintf (stderr, "Line: %d invalid type value: %d\n", __LINE__, type);
    }

    //ret = PyLong_FromLong (rc);
    ret = Py_BuildValue ("ii", rc, new_time);

    return ret;
}

//int stitch_read_block (const StitchFile * file, double timestep, int32_t * dims, int32_t * buffer);
static PyObject * stitch_read_block_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    double time = 0;
    int32_t x1 = 0;
    int32_t y1 = 0;
    int32_t z1 = 0;
    int32_t x2 = 0;
    int32_t y2 = 0;
    int32_t z2 = 0;
    PyObject * dims_parse = NULL;
    int64_t field_id = 0;

    int32_t * passed_dims = NULL;
    int32_t * state_i32 = NULL;
    int64_t * state_i64 = NULL;
    double * state_f64 = NULL;

    int rc = 0;
    PyObject * ret = NULL;
    PyObject * state_ret = NULL;
    int32_t new_time = false;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LLdO", &file, &field_id, &time, &dims_parse))
    {
        return NULL;
    }

    PyArrayObject * dim_obj = (PyArrayObject *) PyArray_FROM_OTF (dims_parse ,NPY_INT32, NPY_ARRAY_FARRAY);
    assert (dim_obj);
    passed_dims = (int32_t *) PyArray_DATA (dim_obj);
    assert (passed_dims);

    x1 = passed_dims [0];
    x2 = passed_dims [1];
    y1 = passed_dims [2];
    y2 = passed_dims [3];
    z1 = passed_dims [4];
    z2 = passed_dims [5];

    enum STITCH_TYPES type = STITCH_NO_TYPE;
    rc = stitch_get_field_type ((StitchFile *) file, field_id, &type);

    // Build py array from scratch with fortran layout
    int32_t nd = 3;
    npy_intp nyp_new_dims [3] = {(x2 - x1), (y2 - y1), (z2 - z1)};
    int flags = NPY_ARRAY_FARRAY;
    PyArray_Descr * state_descr = NULL;

    switch (type)
    {
        case STITCH_INT32:
            state_descr = PyArray_DescrFromType (NPY_INT32);
            state_ret = PyArray_NewFromDescr (&PyArray_Type, state_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
            state_i32 = (int32_t *) PyArray_DATA ((PyArrayObject *) state_ret);
            assert (state);
            break;

        case STITCH_INT64:
            state_descr = PyArray_DescrFromType (NPY_INT64);
            state_ret = PyArray_NewFromDescr (&PyArray_Type, state_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
            state_i64 = (int64_t *) PyArray_DATA ((PyArrayObject *) state_ret);
            assert (state);
            break;

        case STITCH_FLOAT64:
            state_descr = PyArray_DescrFromType (NPY_FLOAT64);
            state_ret = PyArray_NewFromDescr (&PyArray_Type, state_descr, nd, nyp_new_dims, NULL, NULL, flags, NULL);
            state_f64 = (double *) PyArray_DATA ((PyArrayObject *) state_ret);
            assert (state);
            break;

        case STITCH_NO_TYPE:
        default:
            fprintf (stderr, "Line: %d Invalid type: %d\n", __LINE__, type);
    }

    switch (type)
    {
        case STITCH_INT32:
            stitch_read_block_int32 ((StitchFile *) file, field_id, &time, passed_dims, state_i32, &new_time);
            break;

        case STITCH_INT64:
            stitch_read_block_int64 ((StitchFile *) file, field_id, &time, passed_dims, state_i64, &new_time);
            break;

        case STITCH_FLOAT64:
            stitch_read_block_float64 ((StitchFile *) file, field_id, &time, passed_dims, state_f64, &new_time);
            break;

        case STITCH_NO_TYPE:
        default:
            fprintf (stderr, "Line: %d Invalid type: %d\n", __LINE__, type);
    }

#if 0
int32_t min_val = 2000000000;
int32_t max_val = -100;
for (int i = 0; i < ((x2 - x1) * (y2 - y1) * (z2 -z1)); i++)
{
    min_val = (state [i] < min_val ? state [i] : min_val);
    max_val = (state [i] > max_val ? state [i] : max_val);
}

printf ("in stitchmodule: elements: %d min: %d max: %d\n", ((x2 - x1) * (y2 - y1) * (z2 -z1)), min_val, max_val);
#endif

    ret = Py_BuildValue ("iNi", rc, state_ret, new_time);
    assert (ret);

    return ret;
}

// int stitch_get_global_bounds(const StitchFile * file, int64_t field_id, int32_t * bb);
static PyObject * stitch_get_global_bounds_wrapper (PyObject * self, PyObject * args)
{
    int64_t file = 0;
    int64_t field_id = 0;
    // Needs to be 32bit int
    int32_t *block_i32 = NULL;

    int rc = 0;
    PyObject * ret = NULL;
    PyObject * block = NULL;
    PyArray_Descr * descr = NULL;
    // returning fortran layout array
    int flags = NPY_ARRAY_F_CONTIGUOUS;

    // parse arguments
    if (!PyArg_ParseTuple (args, "LL", &file, &field_id))
    {
        return NULL;
    }

    // Build py array from scratch with fortran layout
    // shape=(2,3)
    int32_t nd = 2;
    npy_intp dims [2] = {2,3};
    descr = PyArray_DescrFromType (NPY_INT32);
    block = PyArray_NewFromDescr (&PyArray_Type, descr, nd, dims, NULL, NULL, flags, NULL);
    block_i32 = (int32_t *) PyArray_DATA ((PyArrayObject *) block);

    // Read global bounds
    rc=stitch_get_global_bounds((StitchFile *) file, field_id, block_i32);

    ret = Py_BuildValue ("iN", rc, block);
    assert (ret);

    return ret;
}



//=============================================================================
// Method definition object for this extension, these argumens mean:
// ml_name: The name of the method
// ml_meth: Function pointer to the method implementation
// ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
// ml_doc:  Contents of this method's docstring
static PyMethodDef stitch_methods[] = { 
    {
        "open", stitch_open_wrapper, METH_VARARGS,
        "Open a stitch file for reading and writing"
    },
    {
        "close", stitch_close_wrapper, METH_VARARGS,
        "close a stitch file"
    },
    {
        "set_parameters", stitch_set_parameters_wrapper, METH_VARARGS,
        "set the global parameters for a stitch file"
    },
    {
        "get_parameters", stitch_get_parameters_wrapper, METH_VARARGS,
        "get the global parameters for a stitch file"
    },
    {
        "get_times", stitch_get_times_wrapper, METH_VARARGS,
        "get the list of times written to"
    },
    {
        "get_fields", stitch_get_fields_wrapper, METH_VARARGS,
        "get the list of fields in the file"
    },
    {
        "create_field", stitch_create_field_wrapper, METH_VARARGS,
        "create a new field"
    },
    {
        "set_field_no_value_present", stitch_set_field_no_value_present_wrapper, METH_VARARGS,
        "change the no_value_present value for this field from the default"
    },
    {
        "query_field", stitch_query_field_wrapper, METH_VARARGS,
        "get the field_id for a given label"
    },
    {
        "write_block", stitch_write_block_wrapper, METH_VARARGS,
        "write a block"
    },
    {
        "read_block", stitch_read_block_wrapper, METH_VARARGS,
        "read a block"
    },
    {
        "get_global_bounds", stitch_get_global_bounds_wrapper, METH_VARARGS,
        "read a block"
    },
    {NULL, NULL, 0, NULL}
};

// Module definition
// The arguments of this structure tell Python what to call your extension,
// what it's methods are and where to look for it's method definitions
static struct PyModuleDef stitch_definition = { 
    PyModuleDef_HEAD_INIT,
    "libstitch",
    "A Python3 interface for Stitch",
    -1, 
    stitch_methods
};


// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_libstitch (void)
{
    Py_Initialize();

    // something to do with numpy
    import_array ();
#ifdef STITCH_PARALLEL
    import_mpi4py ();
#endif

    return PyModule_Create(&stitch_definition);
}
