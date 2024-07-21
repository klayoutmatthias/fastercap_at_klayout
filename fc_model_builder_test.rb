
# Run this example with 
#
# ./run_klayout.sh -b -r fc_model_builder_test.rb

require "logger"

#require_relative "fc_model_builder"
load File.join(File.dirname(__FILE__), "fc_model_builder.rb")

logger = Logger.new(STDOUT)
logger.level = Logger::INFO

ly = RBA::Layout::new
net1 = ly.create_cell("Net1")
net2 = ly.create_cell("Net2")
net_cells = [net1, net2]

top = ly.create_cell("TOP")
net_cells.each do |cell|
  top.insert(RBA::CellInstArray::new(cell.cell_index, RBA::Trans::new))
end

m1 = ly.layer(1, 0)  # metal1
v1 = ly.layer(2, 0)  # via1 
m2 = ly.layer(3, 0)  # metal3

net1.shapes(m1).insert(RBA::Box::new(0, 0, 2000, 500))
net1.shapes(v1).insert(RBA::Box::new(1600, 100, 1900, 400))
net1.shapes(m2).insert(RBA::Box::new(1500, 0, 3000, 500))

net2.shapes(m1).insert(RBA::Box::new(-1000, 0, -600, 400))
net2.shapes(v1).insert(RBA::Box::new(-900, 100, -700, 300))
net2.shapes(m2).insert(RBA::Box::new(-1000, 0, -600, 400))

bbox = top.bbox.enlarged(3000)

# generate some nitride liner

rm1_nitride_l = RBA::Region::new(bbox)
rm2_nitride_l = RBA::Region::new(bbox)

fcm = FCModelBuilder::new(3.5, ly.dbu, amax: 0.5, b: 0.5, logger: logger)

fcm.add_material("nit", 4.5)

z = 0.0
h = 0.5
hnit = 0.1

layer = m1
net_cells.each do |cell|
  nn = cell.name
  r = RBA::Region::new(cell.begin_shapes_rec(layer))
  fcm.add_conductor(r, nn, z: z, h: h)
  rnit = r.sized(100)
  fcm.add_dielectric(rnit, "nit", z: z, h: hnit + h)
end

fcm.add_dielectric(rm1_nitride_l, "nit", z: z, h: hnit)

z += h

layer = v1
net_cells.each do |cell|
  nn = cell.name
  r = RBA::Region::new(cell.begin_shapes_rec(layer))
  fcm.add_conductor(r, nn, z: z, h: h)
end

z += h

layer = m2
net_cells.each do |cell|
  nn = cell.name
  r = RBA::Region::new(cell.begin_shapes_rec(layer))
  fcm.add_conductor(r, nn, z: z, h: h)
  rnit = r.sized(100)
  fcm.add_dielectric(rnit, "nit", z: z, h: hnit + h)
end

fcm.add_dielectric(rm2_nitride_l, "nit", z: z, h: hnit)

gen = fcm.generate

# self-check
gen.check

# dump as STL file
gen.dump_stl

# write fastcap file
gen.write_fastcap("test")

