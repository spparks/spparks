
=============
Release Notes
=============

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
