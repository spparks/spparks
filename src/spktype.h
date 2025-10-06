/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

// define integer data types used by SPPARKS and associated size limits

// smallint = variables for on-procesor system (nlocal, nghost, etc)
// tagint = variables for site IDs (nglobal, id[], etc)
// bigint = variables for total system tallies (naccept, nreject, etc)

// smallint must be an int, as defined by C compiler
// tagint can be 32-bit or 64-bit int, must be >= smallint
// bigint can be 32-bit or 64-bit int, must be >= tagint

// MPI_SPK_TAGINT = MPI data type corresponding to a tagint
// MPI_SPK_BIGINT = MPI data type corresponding to a bigint

// NOTE: if your machine/MPI does not support "long long" ints,
//       but only "long" ints, then you will need to change
//       MPI_LONG_LONG to MPI_LONG, and atoll to atol

#ifndef SPK_SPKTYPE_H
#define SPK_SPKTYPE_H

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include "limits.h"
#include "stdint.h"
#include "inttypes.h"

namespace SPPARKS_NS {

// default to 32-bit smallint and tagint, 64-bit bigint

#if !defined(SPPARKS_SMALLSMALL) && !defined(SPPARKS_BIGBIG) && !defined(SPPARKS_SMALLBIG)
#define SPPARKS_SMALLBIG
#endif

// default settings, models restricted to less than 2B sites
// 32-bit smallint and tagint, 64-bit bigint

#ifdef SPPARKS_SMALLBIG

typedef int smallint;
typedef int tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT64_MAX

#define MPI_SPK_TAGINT MPI_INT
#define MPI_SPK_BIGINT MPI_LONG_LONG

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoi
#define ATOBIGINT atoll

#endif

// for models with over 2B sites
// 32-bit smallint, 64-bit tagint and bigint

#ifdef SPPARKS_BIGBIG

typedef int smallint;
typedef int64_t tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT64_MAX
#define MAXBIGINT INT64_MAX

#define MPI_SPK_TAGINT MPI_LONG_LONG
#define MPI_SPK_BIGINT MPI_LONG_LONG

#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoll
#define ATOBIGINT atoll

#endif

// for machines that do not support 64-bit ints
// 32-bit smallint and tagint and bigint

#ifdef SPPARKS_SMALLSMALL

typedef int smallint;
typedef int tagint;
typedef int bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT_MAX

#define MPI_SPK_TAGINT MPI_INT
#define MPI_SPK_BIGINT MPI_INT

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"

#define ATOTAGINT atoi
#define ATOBIGINT atoi

#endif

}

#endif
