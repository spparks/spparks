# viz grain pattern

m = mdump("tmp.grain")
#m = mdump("tmp.membrane")
m.map(2,"spin")
m.etype = "spin"

#e = ensight(m)
#e.one("spin","Spin")

g = gl(m)
v = vcr(g)

