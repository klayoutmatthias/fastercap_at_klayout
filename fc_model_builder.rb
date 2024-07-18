
# @begin A class providing a service for building FastCap2 or FasterCap models
#
# This class is use the following way:
# 
# 1) Create a FCModelBuilder object
#    Specify the default k value which is the k assumed for "empty space".
#    You can also specify a maximum area and the "b" parameter for the 
#    triangulation. The b parameter corresponds to the minimum angle
#    and should be >=1 (b=0.5/sin(a_min)).
#    I.e. b=1 -> a_min=30 deg, b=2 -> a_min~14.5 deg.
#
# 2) Add material definitions for the dielectrics
# 
#    Each material definition consists of a k value and 
#    a material name.
#
# 3) Add layers in the 2.5d view fashion
#    Each layer is a sheet in 3d space that is extruded in vertical 
#    direction with the given start and stop z (or height)
#    The layer must be a DRC::Layer or RBA::Region object.
#    
#    Layers can be added in two ways:
# 
#    * As conductors: specify the net name
#
#    * As dielectric layer: specify the material name
#
# 4) Generate a 3d model using "generate_model"
#    The model is a string with the geometry definitions. This string
#    can be used as input for FasterCap or FastCap2.

# TODO:
# @@@ Generator for the files
# @@@ Some debugging features ..

class FCModelBuilder

  def initialize(k_space, dbu, **kwargs)
    @k_space = k_space
    @dbu = dbu
    @b_value = kwargs[:b] || 1.0
    @amin_value = kwargs[:amin]
    @materials = {}
    @dlayers = {}
    @clayers = {}
    @z = 0.0
    @res = 0.001 # 1nm resolution
  end
  
  def add_material(name, k)
    @materials[name] = k
  end
  
  def add_dielectric(l, mat_name, **kwargs)
    _add_layer(l, true, mat_name, **kwargs)
  end
  
  def add_conductor(l, net_name, **kwargs)
    _add_layer(l, false, net_name, **kwargs)
  end
  
  def _z2norm
    (z / @res + 1e-6).floor
  end
  
  def _norm2z
    z * @res
  end
  
  def _add_layer(l, dielectric, name, **kwargs)
    if dielectric && ! @materials[name]
      raise("Unknown material #{name} - did you use 'add_material'?")
    end
    zstart = kwargs[:z] || kwargs[:zstart] || @z
    if kwargs[:h] 
      zstop = zstart + kwargs[:h]
    else 
      zstop = kwargs[:zz] || kwargs[:zstop] || raise("Stop z has to be given either as height (h) or explicitly (zz/zstop)")
    end
    if dielectric
      (@dlayers[name] ||= []) << [ l, z2norm(zstart), z2norm(zstop) ]
    else
      (@clayers[name] ||= []) << [ l, z2norm(zstart), z2norm(zstop) ]
    end
  end
  
  def generate
  
    z = []
    @layers.each do |k,v|
      v.each do |l|
        z += l[2..3]
      end
    end
    z = z.sort.uniq
    
    if z.empty?
      return
    end

    gen = FCModelGenerator::new(@materials, @clayers.keys)
    gen.k_default = @k_space
    gen.amin = @amin_value
    gen.b = @b_value
    gen.dbu = @dbu
    
    z.each do |zcurr|
    
      gen.next_z(_norm2z(zcurr))
    
      @dlayers.each do |mn,v|
        v.each do |l|
          if l[2] <= zcurr && l[3] > zcurr
            gen.add_in(l[0], "-" + mn)
          end
          if l[2] < zcurr && l[3] >= zcurr
            gen.add_out(l[0], "-" + mn)
          end
        end
      end
      
      @clayers.each do |nn,v|
        v.each do |l|
          if l[2] <= zcurr && l[3] > zcurr
            gen.add_in(l[0], "+" + nn)
          end
          if l[2] < zcurr && l[3] >= zcurr
            gen.add_out(l[0], "+" + nn)
          end
        end
      end
      
      gen.finish_z
    
    end
  
  end

end

class FCModelGenerator

  attr_accessor :k_default
  attr_accessor :amin
  attr_accessor :b
  attr_accessor :dbu

  def initialize(materials, net_names)
    @materials = materials
    @net_names = net_names
    @z = nil
    @zz = nil
    @state = {}
    @k_default = 1.0
    @amin = 0.0
    @b = 2.0
    @diel_data = {}
    @cond_data = {}
    @dbu = 0.001
    self.reset
  end
  
  def reset
    @all_in = RBA::Region::new
    @all_out = RBA::Region::new
    @layers_in = {}
    @layers_out = {}
  end
  
  def add_in(l, name)
    (@layers_in[name] ||= []) << l - @all_in
    @all_in += l
  end
  
  def add_out(l, name)
    (@layers_out[name] ||= []) << l - @all_out
    @all_out += l
  end
  
  def finish_z
  
    din = {}
    dout = {}
    
    all_cin = RBA::Region::new
    all_cout = RBA::Region::new
    @net_names.each do |nn|
      mk = "+" + nn
      lin = @layers_in[mk] || RBA::Region::new
      lout = @layers_in[mk] || RBA::Region::new
      # cancel in/out and out/in
      lin, lout = lin - lout, lout - lin
      all_cin += lin
      all_cout += lout
      @state[mk] += lin
      @state[mk] -= lout
    end
    
    @materials.keys.each do |mn|
      mk = "-" + mn
      lin = @layers_in[mk] || RBA::Region::new
      lout = @layers_in[mk] || RBA::Region::new
      # cancel in/out and out/in
      lin, lout = lin - lout, lout - lin
      lin -= all_cin
      lout = all_cout
      din[mn] = lin
      dout[mn] = lout
      @state[mk] ||= RBA::Region::new
      @state[mk] += lin
      @state[mk] -= lout
    end
    
    @materials.keys.each do |mni|
      lin = din[mni]
      if lin
        lin = lin.dup
        lin -= all_cout  # handled with the conductor
        @materials.keys.each do |mno|
          lout = dout[mno]
          if lout
            d = lout & lin
            if !d.empty?
              generate_hdiel(mno, mni, d)
            end
            lin -= lout
          end
        end
        if !lin.empty?
          generate_hdiel(nil, mni, lin)
        end
      end
    end
    
    @materials.keys.each do |mno|
      lout = dout[mno]
      if lout
        lout = lout.dup
        lout -= all_cin  # handled with the conductor
        @materials.keys.each do |mni|
          lin = din[mni]
          if lin
            lout -= lin
          end
        end
        if !lout.empty?
          generate_hdiel(mno, nil, lout)
        end
      end
    end
    
    @net_names.each do |nn|
      mk = "+" + nn
      lin = @layers_in[mk]
      if lin
        lin = lin.dup
        @materials.keys.each do |mno|
          lout = dout[mno]
          if lout
            d = lout & lin
            if !d.empty?
              generate_hcond_in(nn, mno, d)
            end
            lin -= lout
          end
        end
        if !lin.empty?
          generate_hcond_in(nn, nil, lin)
        end
      end
    end
    
    @net_names.each do |nn|
      mk = "+" + nn
      lout = @layers_out[mk]
      if lout
        lout = lout.dup
        lout -= all_cin  # handled with the conductor
        @materials.keys.each do |mni|
          lin = din[mni]
          if lin
            d = lout & lin
            if !d.empty?
              generate_hcond_out(nn, mni, d)
            end
            lout -= lin
          end
        end
        if !lout.empty?
          generate_hcond_out(nn, nil, lout)
        end
      end
    end
    
  end
  
  def next_z(z)
  
    if !@z
      @z = z
      return
    end
    
    @zz = z
    
    all_cond = RBA::Region::new
    @net_names.each do |nn|
      mk = "+" + nn
      if @state[mk]
        all_cond += @state[mk]
      end
    end
    all_cond = all_cond.edges
    
    @materials.keys.each do |mni|
      linside = @state["-" + mni]
      if linside
        linside = linside.edges
        linside -= all_cond  # handled with the conductor
        @materials.keys.each do |mno|
          if mno != mni
            loutside = @state["-" + mno]
            if loutside
              d = loutside & linside
              d.each do |e|
                generate_vdiel(mno, mni, e)
              end
              linside -= loutside
            end
          end
        end
        linside.each do |e|
          generate_vdiel(nil, mni, e)
        end
      end
    end
    
    @materials.keys.each do |mno|
      loutside = @state["-" + mno]
      if loutside
        loutside = loutside.edges
        loutside -= all_cond  # handled with the conductor
        @materials.keys.each do |mni|
          if mno != mni
            linside = state["-" + mni]
            if linside
              loutside -= linside
            end
          end
        end
        loutside.each do |e|
          generate_vdiel(nn, mno, nil, e)
        end
      end
    end
    
    @net_names.each do |nn|
      mk = "+" + nn
      linside = @state[mk]
      if linside
        linside = linside.edges
        @materials.keys.each do |mno|
          loutside = @state["-" + mno]
          if loutside
            loutside = loutside.edges
            d = loutside & linside
            d.each do |e|
              generate_vcond(nn, mno, e)
            end
            linside -= loutside
          end
        end
        linside.each do |e|
          generate_vcond(nn, nil, e)
        end
      end
    end
    
    @z = z
    
  end  
  
  def generate_hdiel(below, above, layer)
    k = [ below, above ]
    data = (@diel_data[k] ||= [])
    rp = nil
    layer.delaunay(@amin, @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
    end
  end

  def generate_vdiel(left, right, edge)
    el = edge.length
    if el == 0
      return
    end
    k = [ left, right ]
    dbu_trans = RBA::CplxTrans::new(@dbu)
    data = (@diel_data[k] ||= [])
    r = RBA::Region::new
    r.insert(RBA::Box::new(0, 0, el, (@zz - @z) / @dbu))
    layer.delaunay(@amin, @b).each do |t|
      p0 = dbu_trans * edge.p1
      d = dbu_trans * edge.d
      # note: normal is facing outwards (to "left")
      tri = t.each.collect do |pt| 
        pxy = p0 + d * (pt.x.to_f / el.to_f)
        pz = pt.y * @dbu + @z
        [pxy.x, pxy.y, pz]
      end
      data << tri
    end
  end

  def generate_hcond_in(nn, below, layer)
    k = [ nn, below ]
    data = (@cond_data[k] ||= [])
    rp = nil
    layer.delaunay(@amin, @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
    end
  end

  def generate_hcond_out(nn, above, layer)
    k = [ nn, above ]
    data = (@cond_data[k] ||= [])
    rp = nil
    layer.delaunay(@amin, @b).each do |t|
      # note: normal is facing downwards (into conductor)
      tri = t.each.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      # now it is facing outside (to "above")
      data << tri.reverse
    end
  end

  def generate_vcond(nn, left, edge)
    el = edge.length
    if el == 0
      return
    end
    k = [ nn, left ]
    dbu_trans = RBA::CplxTrans::new(@dbu)
    data = (@cond_data[k] ||= [])
    r = RBA::Region::new
    r.insert(RBA::Box::new(0, 0, el, (@zz - @z) / @dbu))
    layer.delaunay(@amin, @b).each do |t|
      p0 = dbu_trans * edge.p1
      d = dbu_trans * edge.d
      # note: normal is facing outwards (to "left")
      tri = t.each.collect do |pt| 
        pxy = p0 + d * (pt.x.to_f / el.to_f)
        pz = pt.y * @dbu + @z
        [pxy.x, pxy.y, pz]
      end
      data << tri
    end
  end

end

