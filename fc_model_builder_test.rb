
#require_relative "fc_model_builder"
load File.join(File.dirname(__FILE__), "fc_model_builder.rb")

ly = RBA::Layout::new
top = ly.create_cell("TOP")

m1 = ly.layer(1, 0)  # metal1
v1 = ly.layer(2, 0)  # via1 
m2 = ly.layer(3, 0)  # metal3

top.shapes(m1).insert(RBA::Box::new(0, 0, 2000, 500))
top.shapes(v1).insert(RBA::Box::new(1600, 100, 1900, 400))
top.shapes(m2).insert(RBA::Box::new(1500, 0, 3000, 500))

bbox = top.bbox.enlarged(3000)

rm1 = RBA::Region::new(top.begin_shapes_rec(m1))
rv1 = RBA::Region::new(top.begin_shapes_rec(v1))
rm2 = RBA::Region::new(top.begin_shapes_rec(m2))

# generate some nitride liner

rm1_nitride_l = RBA::Region::new(bbox)
rm1_nitride_h = rm1.sized(100)

rm2_nitride_l = RBA::Region::new(bbox)
rm2_nitride_h = rm2.sized(100)

fcm = FCModelBuilder::new(3.5, ly.dbu, amin: 0.1, b: 0.5)

fcm.add_material("nit", 7.0)

z = 0.0
h = 0.5
hnit = 0.1
fcm.add_conductor(rm1, "Net1", z: z, h: h)
fcm.add_dielectric(rm1_nitride_l, "nit", z: z, h: hnit)
fcm.add_dielectric(rm1_nitride_h, "nit", z: z, h: hnit + h)

z += h
fcm.add_conductor(rv1, "Net1", z: z, h: h)

z += h
fcm.add_conductor(rm2, "Net1", z: z, h: h)
fcm.add_dielectric(rm2_nitride_l, "nit", z: z, h: hnit)
fcm.add_dielectric(rm2_nitride_h, "nit", z: z, h: hnit + h)

fcm.generate
