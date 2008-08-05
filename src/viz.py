# viz grain pattern

d = dump("tmp.pore")
d.aselect.test("$type == 2")
g = gl(d)
g.acol(2,"blue")
g.arad(0,0.5)
v = vcr(g)
