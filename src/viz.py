# viz grain pattern

d = dump("tmp.dep")
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")
d.aselect.test("$type == 2")
g = gl(d)
g.acol(2,"blue")
v = vcr(g)

#d = dump("tmp.membrane")
#g = gl(d)
#g.arad(0,0.5)
#v = vcr(g)
