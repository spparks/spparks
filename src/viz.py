# viz grain pattern

d = dump("tmp.pore2")
d.aselect.test("$type == 2")
g = gl(d)
g.acol(2,"blue")
v = vcr(g)

#d = dump("tmp.membrane")
#g = gl(d)
#g.arad(0,0.5)
#v = vcr(g)
