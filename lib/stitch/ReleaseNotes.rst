
=============
Release Notes
=============

May 2022
--------

This stitch release includes a new 'global-bounds' functionality.  This is
additional capability allows user to query and existing stitch database for a
particular field: 'what are the global bounds' written to the file for this
field across all time steps; a bounding box is returned; this can be useful 
when workflows use a stitch file that was written upstream and knowledge of 
the spatial extent of the file is not known or it is inconvenient to carry 
that information.  Note that this global bounds is not stored in the stitch 
file as meta data -- rather it is computed on the fly.

Jan 18 2021
-----------

This is a v1.1 stitch update which substantially improves performance of the
'set stitch' command in spparks. The 'set stitch' command is most commonly used
for additive manufacturing (AM) simulations where many spparks runs are
'stitched' together; each run requires use of the 'set stitch' command; as the
AM simulation proceeds, stitch files grows proportionally in size which
increases IO time relating to stitch database queries. This release improves
the performance of stitch queries; query times now scale linearly with stitch
file size.  There are no new files with this release.  The algorithm used in
stitch query has been substantially improved.  For small problems, on the order
of millions of lattice sites, the performance improvements may not be
noticable.  For problem sizes on the order 10s to 100s of millions of lattice
sites performance gains are significant.
