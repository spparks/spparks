# viz grain pattern

#d = dump("tmp.pore")
#d.aselect.test("$type == 2")
#g.acol(2,"blue")

d = dump("tmp.membrane")
g = gl(d)
g.arad(0,0.5)
v = vcr(g)
