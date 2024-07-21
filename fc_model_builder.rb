
# @begin A class providing a service for building FastCap2 or FasterCap models
#
# This class is used the following way:
# 
# 1) Create a FCModelBuilder object
#    Specify the default k value which is the k assumed for "empty space".
#    You can also specify a maximum area and the "b" parameter for the 
#    triangulation. The b parameter corresponds to the minimum angle
#    and should be <=1 (b=sin(min_angle)*2).
#    I.e. b=1 -> min_angle=30 deg, b=0.5 -> min_angle~14.5 deg.
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
    @amax_value = kwargs[:amax] || 0.0
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
    if l.respond_to?(:data)
      l = l.data
    end
    _add_layer(l, true, mat_name, **kwargs)
  end
  
  def add_conductor(l, net_name, **kwargs)
    if l.respond_to?(:data)
      l = l.data
    end
    _add_layer(l, false, net_name, **kwargs)
  end
  
  def _z2norm(z)
    (z / @res + 1e-6).floor
  end
  
  def _norm2z(z)
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
      (@dlayers[name] ||= []) << [ l, _z2norm(zstart), _z2norm(zstop) ]
    else
      (@clayers[name] ||= []) << [ l, _z2norm(zstart), _z2norm(zstop) ]
    end
  end
  
  def generate
  
    z = []
    [ @dlayers, @clayers ].each do |ll|
      ll.each do |k,v|
        v.each do |l|
          z += l[1..2]
        end
      end
    end
    z = z.sort.uniq
    
    if z.empty?
      return
    end

    gen = FCModelGenerator::new(@materials, @clayers.keys)
    gen.k_default = @k_space
    gen.amax = @amax_value
    gen.b = @b_value
    gen.dbu = @dbu
    
    z.each do |zcurr|
    
      gen.next_z(_norm2z(zcurr))
    
      @clayers.each do |nn,v|
        v.each do |l|
          if l[1] <= zcurr && l[2] > zcurr
            gen.add_in(l[0], "+" + nn)
          end
          if l[1] < zcurr && l[2] >= zcurr
            gen.add_out(l[0], "+" + nn)
          end
        end
      end
      
      @dlayers.each do |mn,v|
        v.each do |l|
          if l[1] <= zcurr && l[2] > zcurr
            gen.add_in(l[0], "-" + mn)
          end
          if l[1] < zcurr && l[2] >= zcurr
            gen.add_out(l[0], "-" + mn)
          end
        end
      end
      
      gen.finish_z
    
    end

    # self-check
    gen.check
  
  end

end

class FCModelGenerator

  attr_accessor :k_default
  attr_accessor :amax
  attr_accessor :b
  attr_accessor :dbu

  def initialize(materials, net_names)
    @materials = materials
    @net_names = net_names
    @z = nil
    @zz = nil
    @state = {}
    @k_default = 1.0
    @amax = 0.0
    @b = 0.5
    @diel_data = {}
    @cond_data = {}
    @dbu = 0.001
  end
  
  def reset
    @layers_in = {}
    @layers_out = {}
  end
  
  def add_in(l, name)
    puts "add_in: #{name} -> #{l.to_s}" 
    @layers_in[name] ||= RBA::Region::new
    @layers_in[name] += l
  end
  
  def add_out(l, name)
    puts "add_out: #{name} -> #{l.to_s}" 
    @layers_out[name] ||= RBA::Region::new
    @layers_out[name] += l
  end
  
  def finish_z
  
    puts "Finishing layer z=#{@z}"

    din = {}
    dout = {}
    
    all_in = RBA::Region::new
    all_out = RBA::Region::new
    all = RBA::Region::new

    all_cin = nil
    all_cout = nil

    [ [ @net_names, "+" ], [ @materials.keys, "-" ] ].each do |np|

      names, prefix = np

      names.each do |nn|

        mk = prefix + nn
        lin = @layers_in[mk] ? @layers_in[mk].dup : RBA::Region::new
        lout = @layers_out[mk] ? @layers_out[mk].dup : RBA::Region::new

        # cancel in/out and out/in
        lin, lout = lin - lout, lout - lin

        @state[mk] ||= RBA::Region::new
        st = @state[mk] + lin - lout

        lout -= all
        lin -= all
        lout += st & all_in
        lin += st & all_out

        if @layers_out[mk]
          lin -= @layers_out[mk]
        end
        if @layers_in[mk]
          lout -= @layers_in[mk]
        end

        @state[mk] += lin
        @state[mk] -= lout

        din[mk] = lin
        dout[mk] = lout

        puts "Legalized events & status for '#{mk}':"
        puts "  in = #{din[mk].to_s}"
        puts "  out = #{dout[mk].to_s}"
        puts "  state = #{@state[mk].to_s}"

        all_in += lin
        all_out += lout
        all += @state[mk]

      end

      if prefix == "+"
        all_cin = all_in.dup
        all_cout = all_out.dup
      end

    end

    puts "All conductor region in: #{all_cin.to_s}"
    puts "All conductor region out: #{all_cout.to_s}"

    # check whether the states are separated
    a = @state.values.inject(:+)
    @state.each do |k,s|
      r = s - a
      if ! r.is_empty?
        puts "*** ERROR: state region of '#{k}' (#{s.to_s}) is not contained entirely in remaining all state region (#{a.to_s}) - this means there is an overlap"
      end
      a -= s
    end

    # Now we have legalized the in and out events
    
    @materials.keys.each do |mni|
      lin = din["-" + mni]
      if lin
        lin = lin.dup
        lin -= all_cout  # handled with the conductor
        @materials.keys.each do |mno|
          lout = dout["-" + mno]
          if lout
            d = lout & lin
            if !d.is_empty?
              generate_hdiel(mno, mni, d)
            end
            lin -= lout
          end
        end
        if !lin.is_empty?
          generate_hdiel(nil, mni, lin)
        end
      end
    end
    
    @materials.keys.each do |mno|
      lout = dout["-" + mno]
      if lout
        lout = lout.dup
        lout -= all_cin  # handled with the conductor
        @materials.keys.each do |mni|
          lin = din["-" + mni]
          if lin
            lout -= lin
          end
        end
        if !lout.is_empty?
          generate_hdiel(mno, nil, lout)
        end
      end
    end
    
    @net_names.each do |nn|
      lin = din["+" + nn]
      if lin
        lin = lin.dup
        @materials.keys.each do |mno|
          lout = dout["-" + mno]
          if lout
            d = lout & lin
            if !d.is_empty?
              generate_hcond_in(nn, mno, d)
            end
            lin -= lout
          end
        end
        if !lin.is_empty?
          generate_hcond_in(nn, nil, lin)
        end
      end
    end
    
    @net_names.each do |nn|
      lout = dout["+" + nn]
      if lout
        lout = lout.dup
        lout -= all_cin  # handled with the conductor
        @materials.keys.each do |mni|
          lin = din["-" + mni]
          if lin
            d = lout & lin
            if !d.is_empty?
              generate_hcond_out(nn, mni, d)
            end
            lout -= lin
          end
        end
        if !lout.is_empty?
          generate_hcond_out(nn, nil, lout)
        end
      end
    end
    
  end
  
  def next_z(z)
  
    puts "Next layer z=#{z} .."

    self.reset 

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
    puts "Generating horizontal dielectric surface #{below || '(void)'} <-> #{above || '(void)'} as #{layer.to_s}"
    k = [ below, above ]
    data = (@diel_data[k] ||= [])
    rp = nil
    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
    end
  end

  def generate_vdiel(left, right, edge)
    puts "Generating vertical dielectric surface #{left || '(void)'} <-> #{right || '(void)'} with edge #{edge.to_s}"
    el = edge.length
    if el == 0
      return
    end
    k = [ left, right ]
    dbu_trans = RBA::CplxTrans::new(@dbu)
    data = (@diel_data[k] ||= [])
    r = RBA::Region::new
    r.insert(RBA::Box::new(0, 0, el, (@zz - @z) / @dbu))
    r.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      p0 = dbu_trans * edge.p1
      d = dbu_trans * edge.d
      # note: normal is facing outwards (to "left")
      tri = t.each_point_hull.collect do |pt| 
        pxy = p0 + d * (pt.x.to_f / el.to_f)
        pz = pt.y * @dbu + @z
        [pxy.x, pxy.y, pz]
      end
      data << tri
    end
  end

  def generate_hcond_in(nn, below, layer)
    puts "Generating horizontal bottom conductor surface #{below || '(void)'} <-> #{nn} as #{layer.to_s}"
    k = [ nn, below ]
    data = (@cond_data[k] ||= [])
    rp = nil
    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
    end
  end

  def generate_hcond_out(nn, above, layer)
    puts "Generating horizontal top conductor surface #{nn} <-> #{above || '(void)'} as #{layer.to_s}"
    k = [ nn, above ]
    data = (@cond_data[k] ||= [])
    rp = nil
    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (into conductor)
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      # now it is facing outside (to "above")
      data << tri.reverse
    end
  end

  def generate_vcond(nn, left, edge)
    puts "Generating vertical conductor surface #{nn} <-> #{left || '(void)'} with edge #{edge.to_s}"
    el = edge.length
    if el == 0
      return
    end
    k = [ nn, left ]
    dbu_trans = RBA::CplxTrans::new(@dbu)
    data = (@cond_data[k] ||= [])
    r = RBA::Region::new
    r.insert(RBA::Box::new(0, 0, el, (@zz - @z) / @dbu))
    r.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      p0 = dbu_trans * edge.p1
      d = dbu_trans * edge.d
      # note: normal is facing outwards (to "left")
      tri = t.each_point_hull.collect do |pt| 
        pxy = p0 + d * (pt.x.to_f / el.to_f)
        pz = pt.y * @dbu + @z
        [pxy.x, pxy.y, pz]
      end
      data << tri
    end
  end

  def check

    puts "Checking .."
    errors = 0

    @materials.each do |mn|

      tris = _collect_diel_tris(mn)
      puts "Material #{mn} -> #{tris.size} triangles"

      edge_hash = {}
      edges = _normed_edges(tris)

      edges.each do |e|
        if edge_hash[e]
          puts "ERROR: material '#{mn}': duplicate edge #{_edge2s(e)}"
          errors += 1
        else
          edge_hash[e] = true
        end
      end 

      edges.each do |e|
        if !edge_hash[e.reverse]
          puts "ERROR: material '#{mn}': edge #{_edge2s(e)} not connected with reverse edge (open surface)"
          errors += 1
        end
      end

    end

    errors

  end

  def _edge2s(e)
    "(%.12g,%.12g,%12g)-(%12g,%12g,%12g)" % (e.p1 + e.p2)
  end

  def _normed_edges(tris)

    edges = []

    tris.each do |t|
      3.times each do |i|
        p1 = t[i]
        p2 = t[(i + 1) % 3]
        p1 = p1.collect { |c| (c / @dbu + 0.5).floor }
        p2 = p2.collect { |c| (c / @dbu + 0.5).floor }
        edges << [ p1, p2 ]
      end
    end

    edges

  end

  def _collect_diel_tris(mn)

    tris = []

    @diel_data.each do |k,v|
      outside, inside = k
      if outside == mn
        tris += v
      elsif inside == mn
        tris += v.collect { |t| t.reverse }
      end
    end
        
    @cond_data.each do |k,v|
      nn, outside = k
      if outside == mn
        tris += v
      end
    end
        
    return tris

  end

  def _collect_cond_tris(net_name)

    tris = []

    @cond_data.each do |k,v|
      nn, outside = k
      if nn == net_name
        tris += v.collect { |t| t.reverse }
      end
    end
        
    return tris

  end

end

