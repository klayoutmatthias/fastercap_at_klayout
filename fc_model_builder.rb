
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
#    The layers can intersect. The package resolves intersections
#    based on priority: conductors first, dielectrics according to
#    their position in the "materials" definition (first entries have
#    higher prio)
#
# 4) Generate a 3d model using "generate"
#    This method returns an object you can use to generate STL files
#    or FastCap files.

class FCModelBuilder

  def initialize(k_void, dbu, **kwargs)

    @k_void = k_void
    @dbu = dbu
    @b_value = kwargs[:b] || 1.0
    @amax_value = kwargs[:amax] || 0.0
    @materials = {}
    @dlayers = {}
    @clayers = {}
    @z = nil
    @logger = kwargs[:logger]

    @logger && @logger.info("DBU: #{'%.12g' % @dbu}")
    @logger && @logger.info("Delaunay b: #{'%.12g' % @b_value}")
    @logger && @logger.info("Delaunay amax: #{'%.12g' % @amax_value} um^2")

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
    (z / @dbu + 1e-6).floor
  end
  
  def _norm2z(z)
    z * @dbu
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

    gen = FCModelGenerator::new(@materials, @clayers.keys, @logger)
    gen.k_void = @k_void
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

    gen.finalize

    gen
  
  end

end

class FCModelGenerator

  attr_accessor :k_void
  attr_accessor :amax
  attr_accessor :b
  attr_accessor :dbu

  def initialize(materials, net_names, logger)

    @logger = logger
    @materials = materials
    @net_names = net_names
    @z = nil
    @zz = nil
    @state = {}
    @current = {}
    @k_void = 1.0
    @amax = 0.0
    @b = 0.5
    @diel_data = {}
    @diel_vdata = {}
    @cond_data = {}
    @cond_vdata = {}
    @dbu = 0.001

  end
  
  def reset
    @layers_in = {}
    @layers_out = {}
  end
  
  def add_in(l, name)
    @logger && @logger.debug("add_in: #{name} -> #{l.to_s}" )
    @layers_in[name] ||= RBA::Region::new
    @layers_in[name] += l
  end
  
  def add_out(l, name)
    @logger && @logger.debug("add_out: #{name} -> #{l.to_s}" )
    @layers_out[name] ||= RBA::Region::new
    @layers_out[name] += l
  end
  
  def finish_z
  
    @logger && @logger.debug("Finishing layer z=#{@z}")

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

        # compute merged events
        @current[mk] ||= []
        current_before = @current[mk][0] ? @current[mk][0].dup : RBA::Region::new
        lin, lout, current = _merge_events(@current[mk], @layers_in[mk], @layers_out[mk])
        @logger && @logger.debug("Merged events & status for '#{mk}':")
        @logger && @logger.debug("  in = #{lin.to_s}")
        @logger && @logger.debug("  out = #{lout.to_s}")
        @logger && @logger.debug("  state = #{current.to_s}")

        @state[mk] ||= RBA::Region::new
        
        # legalize in and out events
        lin_org = lin.dup
        lout_org = lout.dup
        lout &= @state[mk]
        lin -= all
        lout += current & all_in
        lin += current_before & all_out
        lin -= lout_org
        lout -= lin_org

        # tracks the legalized horizontal cuts
        @state[mk] += lin
        @state[mk] -= lout

        din[mk] = lin
        dout[mk] = lout

        @logger && @logger.debug("Legalized events & status for '#{mk}':")
        @logger && @logger.debug("  in = #{din[mk].to_s}")
        @logger && @logger.debug("  out = #{dout[mk].to_s}")
        @logger && @logger.debug("  state = #{@state[mk].to_s}")

        all_in += lin
        all_out += lout
        all += @state[mk]

      end

      if prefix == "+"
        all_cin = all_in.dup
        all_cout = all_out.dup
      end

    end

    @logger && @logger.debug("All conductor region in: #{all_cin.to_s}")
    @logger && @logger.debug("All conductor region out: #{all_cout.to_s}")

    # check whether the states are separated
    a = @state.values.inject(:+)
    @state.each do |k,s|
      r = s - a
      if ! r.is_empty?
        @logger && @logger.error("State region of '#{k}' (#{s.to_s}) is not contained entirely in remaining all state region (#{a.to_s}) - this means there is an overlap")
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
  
    @logger && @logger.debug("Next layer z=#{z} ..")

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
    
    mat_keys = @materials.keys
    mat_keys.size.times do |i|
      mni = mat_keys[i]
      linside = @state["-" + mni]
      if linside
        linside = linside.edges
        linside -= all_cond  # handled with the conductor
        mat_keys.size.times do |o|
          if i != o
            mno = mat_keys[o]
            loutside = @state["-" + mno]
            if loutside
              loutside = loutside.edges
              if o > i
                d = loutside & linside
                d.each do |e|
                  # NOTE: we need to swap points as we started from "outside"
                  generate_vdiel(mno, mni, e.swapped_points)
                end
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
              # NOTE: we need to swap points as we started from "outside"
              generate_vcond(nn, mno, e.swapped_points)
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

    @logger && @logger.debug("Generating horizontal dielectric surface #{below || '(void)'} <-> #{above || '(void)'} as #{layer.to_s}")

    k = [ below, above ]
    data = (@diel_data[k] ||= [])

    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
      @logger && @logger.debug("  #{tri.to_s}")
    end

  end

  def generate_vdiel(left, right, edge)

    @logger && @logger.debug("Generating vertical dielectric surface #{left || '(void)'} <-> #{right || '(void)'} with edge #{edge.to_s}")

    if edge.is_degenerate?
      return
    end

    el = Math.sqrt(edge.sq_length)

    de = RBA::DVector::new(edge.d.x / el, edge.d.y / el)
    ne = RBA::DVector::new(edge.d.y / el, -edge.d.x / el)
    p0 = ne * ne.sprod(RBA::DPoint::new(edge.p1) - RBA::DPoint::new) + RBA::DPoint::new

    x1 = (edge.p1 - p0).sprod(de)
    x2 = (edge.p2 - p0).sprod(de)

    k = [ left, right, p0, de ]
    @diel_vdata[k] ||= RBA::Region::new
    @diel_vdata[k].insert(RBA::Box::new(x1, (@z / @dbu + 0.5).floor, x2, (@zz / @dbu + 0.5).floor))

  end

  def generate_hcond_in(nn, below, layer)

    @logger && @logger.debug("Generating horizontal bottom conductor surface #{below || '(void)'} <-> #{nn} as #{layer.to_s}")

    k = [ nn, below ]
    data = (@cond_data[k] ||= [])

    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (to "below")
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      data << tri
      @logger && @logger.debug("  #{tri.to_s}")
    end

  end

  def generate_hcond_out(nn, above, layer)

    @logger && @logger.debug("Generating horizontal top conductor surface #{nn} <-> #{above || '(void)'} as #{layer.to_s}")

    k = [ nn, above ]
    data = (@cond_data[k] ||= [])

    layer.delaunay(@amax / (@dbu * @dbu), @b).each do |t|
      # note: normal is facing downwards (into conductor)
      tri = t.each_point_hull.collect { |pt| [pt.x * @dbu, pt.y * @dbu, @z] }
      # now it is facing outside (to "above")
      data << tri.reverse
      @logger && @logger.debug("  #{tri.reverse.to_s}")
    end

  end

  def generate_vcond(nn, left, edge)

    @logger && @logger.debug("Generating vertical conductor surface #{nn} <-> #{left || '(void)'} with edge #{edge.to_s}")

    if edge.is_degenerate?
      return
    end

    el = Math.sqrt(edge.sq_length)

    de = RBA::DVector::new(edge.d.x / el, edge.d.y / el)
    ne = RBA::DVector::new(edge.d.y / el, -edge.d.x / el)
    p0 = ne * ne.sprod(RBA::DPoint::new(edge.p1) - RBA::DPoint::new) + RBA::DPoint::new

    x1 = (edge.p1 - p0).sprod(de)
    x2 = (edge.p2 - p0).sprod(de)

    k = [ nn, left, p0, de ]
    @cond_vdata[k] ||= RBA::Region::new
    @cond_vdata[k].insert(RBA::Box::new(x1, (@z / @dbu + 0.5).floor, x2, (@zz / @dbu + 0.5).floor))

  end

  def finalize

    dbu_trans = RBA::CplxTrans::new(@dbu)

    @diel_vdata.keys.each do |k|

      r = @diel_vdata[k]
      left, right, p0, de = k

      @logger && @logger.debug("Finishing vertical dielectric plane #{left || 'void'} <-> #{right || 'void'} at #{p0.to_s}/#{de.to_s}")

      kk = [ left, right ]
      data = (@diel_data[kk] ||= [])

      r.delaunay(@amax / (@dbu * @dbu), @b).each do |t|

        # note: normal is facing outwards (to "left")
        tri = t.each_point_hull.collect do |pt| 
          pxy = (p0 + de * pt.x) * @dbu
          pz = pt.y * @dbu
          [pxy.x, pxy.y, pz]
        end

        data << tri

        @logger && @logger.debug("  #{tri.to_s}")

      end

    end

    @cond_vdata.keys.each do |k|

      r = @cond_vdata[k]
      nn, left, p0, de = k

      @logger && @logger.debug("Finishing vertical conductor plane #{nn} <-> #{left || 'void'} at #{p0.to_s}/#{de.to_s}")

      kk = [ nn, left ]
      data = (@cond_data[kk] ||= [])

      r.delaunay(@amax / (@dbu * @dbu), @b).each do |t|

        # note: normal is facing outwards (to "left")
        tri = t.each_point_hull.collect do |pt| 
          pxy = (p0 + de * pt.x) * @dbu
          pz = pt.y * @dbu
          [pxy.x, pxy.y, pz]
        end

        data << tri

        @logger && @logger.debug("  #{tri.to_s}")

      end

    end

    dk = {}
    @diel_data.keys.each do |k|
      kk = [k[1], k[0]]
      if !dk[kk]
        dk[k] = []
      else 
        @logger && @logger.debug("Combining dielectric surfaces #{kk[0] || 'void'} <-> #{kk[1] || 'void'} with reverse")
      end
    end

    @diel_data.each do |k,v|
      kk = [k[1], k[0]]
      if dk[kk]
        dk[kk] += v.collect { |t| t.reverse } 
      else
        dk[k] += v
      end
    end
        
    @diel_data = dk

  end

  def dump_stl

    @materials.keys.each do |mn|

      tris = _collect_diel_tris(mn)
      _write_as_stl("diel_#{mn}.stl", tris)

    end

    @net_names.each do |nn|

      tris = _collect_cond_tris(nn)
      _write_as_stl("cond_#{nn}.stl", tris)

    end

  end

  def write_fastcap(prefix)

    lst_fn = prefix + ".lst"

    file_num = 0
    lst_file = []

    lst_file << "* k_void=" + ("%.12g" % @k_void)

    @diel_data.keys.each do |k|

      data = @diel_data[k]

      if data.empty?
        next
      end

      file_num += 1

      outside, inside = k

      k_outside = outside ? @materials[outside] : @k_void
      k_inside = inside ? @materials[inside] : @k_void

      lst_file << "* Dielectric interface: outside=#{outside || '(void)'}, inside=#{inside || '(void)'}"

      fn = prefix + "_" + file_num.to_s + ".geo"
      _write_fastcap_geo(file_num, fn, data, nil)

      lst_file << "D #{fn}  #{'%.12g' % k_outside}  #{'%.12g' % k_inside}  0 0 0"

    end

    prev_net = nil
    prev_cline = nil
    net_num = 0

    @cond_data.keys.sort { |a,b| a[0] <=> b[0] }.each do |k|

      data = @cond_data[k]

      if data.empty?
        next
      end

      nn, outside = k

      if prev_net == nn
        lst_file[prev_cline] += "  +"
      else
        net_num += 1
      end

      file_num += 1

      k_outside = outside ? @materials[outside] : @k_void

      lst_file << "* Conductor interface: outside=#{outside || '(void)'}, net=#{nn}"

      fn = prefix + "_" + file_num.to_s + ".geo"
      _write_fastcap_geo(net_num, fn, data, nn)

      cline = "C #{fn}  #{'%.12g' % k_outside}  0 0 0"
      prev_cline = lst_file.size
      lst_file << cline

      prev_net = nn

    end

    @logger && @logger.info("Writing FasterCap list file: #{lst_fn}")
    File.open(lst_fn, "w") do |file|
      file.write(lst_file.join("\n"))
      file.write("\n")
    end

  end

  def _write_fastcap_geo(num, fn, data, cond_name)

    @logger && @logger.info("Writing FasterCap geo file: #{fn}")
    File.open(fn, "w") do |file|

      data.each do |t|

        file.write("T #{num}")

        t.each do |p|
          file.write("  " + p.collect { |c| '%.12g' % c }.join(" "))
        end

        # compute a reference point in the "outside"
        dx1, dy1, dz1 = 3.times.collect { |i| t[1][i] - t[0][i] }
        dx2, dy2, dz2 = 3.times.collect { |i| t[2][i] - t[0][i] }
        n = [ dy1 * dz2 - dy2 * dz1, -(dx1 * dz2 - dx2 * dz1), dx1 * dy2 - dx2 * dy1 ]
        n = 3.times.collect { |i| n[i] + t[0][i] }

        # write this reference point
        file.write("  " + n.collect { |c| '%.12g' % c }.join(" "))

        file.write("\n")

      end

      if cond_name
        file.write("N #{num} #{cond_name}\n")
      end

    end

  end

  def _check_tris(msg, tris)

    errors = 0
    
    edge_hash = {}
    edges = _normed_edges(tris)

    edges.each do |e|
      if edge_hash[e]
        @logger && @logger.error(msg + ": duplicate edge #{_edge2s(e)}")
        errors += 1
      else
        edge_hash[e] = true
      end
    end
    
    _split_edges(edge_hash)

    edge_hash.keys.each do |e|
      if !edge_hash[e.reverse]
        @logger && @logger.error(msg + ": edge #{_edge2s(e)} not connected with reverse edge (open surface)")
        errors += 1
      end
    end
    
    errors
  
  end
  
  def _normed_edges(tris)

    edges = []

    tris.each do |t|
      3.times do |i|
        p1 = t[i]
        p2 = t[(i + 1) % 3]
        p1 = p1.collect { |c| (c / @dbu + 0.5).floor }
        p2 = p2.collect { |c| (c / @dbu + 0.5).floor }
        edges << [ p1, p2 ]
      end
    end
    
    edges

  end
  
  def _vector_of_edge(e)
    3.times.collect { |i| e[1][i] - e[0][i] }
  end
  
  def _is_antiparallel(a, b)
    vp = 3.times.collect do |i|
      i1 = (i + 1) % 3
      i2 = (i + 2) % 3
      a[i1] * b[i2] - a[i2] * b[i1]
    end
    if _sq_length(vp).abs > 0.5 # we got normalized!
      return false
    end
    sp = 3.times.collect do |i|
      a[i] * b[i]
    end.inject(:+)
    return sp < 0
  end
  
  def _sq_length(e)
    e.collect { |c| c*c }.inject(:+)
  end
  
  # @brief Splits edges into pieces to match an antiparallel sequence of multiple colinear edges
  #
  # Such cases arise at the connection line between two sub-surfaces which are triangulated
  # separately, because their triangles are not neccessarily aligned at the corners.
  
  def _split_edges(edge_hash)
  
    edges_by_p2 = {}
    edges_by_p1 = {}
    edge_hash.keys.each do |e|
      edges_by_p2[e[1]] ||= []
      edges_by_p2[e[1]] << e
      edges_by_p1[e[0]] ||= []
      edges_by_p1[e[0]] << e
    end
    
    while(true)
    
      subst = {}
    
      edge_hash.keys.each do |e|
        ee = edges_by_p2[e[0]] || []
        (edges_by_p2[e[0]] || []).each do |ee|
          if _is_antiparallel(_vector_of_edge(e), _vector_of_edge(ee)) && _sq_length(_vector_of_edge(ee)) < _sq_length(_vector_of_edge(e)) - 0.5
            # There is a shorter edge antiparallel -> this means we need to insert a split point into e
            subst[e] ||= []
            subst[e] += [ [ e[0], ee[0] ], [ ee[0], e[1] ] ]
          end
        end
      end
      
      edge_hash.keys.each do |e|
        ee = edges_by_p1[e[1]] || []
        (edges_by_p1[e[1]] || []).each do |ee|
          if _is_antiparallel(_vector_of_edge(e), _vector_of_edge(ee)) && _sq_length(_vector_of_edge(ee)) < _sq_length(_vector_of_edge(e)) - 0.5
            # There is a shorter edge antiparallel -> this means we need to insert a split point into e
            subst[e] ||= []
            subst[e] += [ [ e[0], ee[1] ], [ ee[1], e[1] ] ]
          end
        end
      end
      
      if subst.empty?
        break
      end
      
      subst.each do |e,with|
        edges_by_p1[e[0]].delete(e)
        edges_by_p2[e[1]].delete(e)
        edge_hash.delete(e)
        with.each do |w|
          edge_hash[w] = true
          (edges_by_p1[w[0]] ||= []) << w
          (edges_by_p2[w[1]] ||= []) << w
        end
      end
    
    end
  
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

  def check

    @logger && @logger.info("Checking ..")
    errors = 0

    @materials.keys.each do |mn|

      tris = _collect_diel_tris(mn)
      @logger && @logger.info("Material #{mn} -> #{tris.size} triangles")
      
      errors += _check_tris("Material '#{mn}'", tris)

    end

    @net_names.each do |nn|

      tris = _collect_cond_tris(nn)
      @logger && @logger.info("Net #{nn} -> #{tris.size} triangles")

      errors += _check_tris("Net '#{nn}'", tris)

    end

    if errors == 0
      @logger && @logger.info("  No errors found!")
    else
      @logger && @logger.info("  #{errors} error(s) found")
    end

    errors

  end

  def _write_as_stl(filename, tris)
  
    if tris.empty?
      return
    end

    @logger && @logger.info("Writing STL file #{filename} ..")

    File.open(filename, "w") do |file|

      file.write("solid stl\n")

      tris.each do |t|

        file.write(" facet normal 0 0 0\n")
        file.write("  outer loop\n")
        t.reverse.each do |p|
          file.write("   vertex %.12g %.12g %.12g\n" % p)
        end
        file.write("  endloop\n")
        file.write(" endfacet\n")

      end

      file.write("endsolid stl\n")

    end

  end

  def _merge_events(pyra, lin, lout)
    
    if lin
      lin = lin.dup
    else
      lin = RBA::Region::new
    end

    if lout
      lout = lout.dup
    else
      lout = RBA::Region::new
    end

    past = (pyra[0] ? pyra[0].dup : RBA::Region::new)

    pyra.size.times do |i|
      ii = pyra.size - i
      added = lin & pyra[ii - 1]
      if !added.is_empty?
        pyra[ii] ||= RBA::Region::new
        pyra[ii] += added
        lin -= added
      end
    end

    pyra[0] ||= RBA::Region::new
    pyra[0] += lin

    pyra.size.times do |i|
      ii = pyra.size - i
      removed = lout & pyra[ii - 1]
      if !removed.is_empty?
        pyra[ii - 1] -= removed
        lout -= removed
      end
    end

    # compute merged events
    lin = pyra[0] - past
    lout = past - pyra[0]

    return [lin, lout, pyra[0]]

  end

  def _edge2s(e)
    "(%.12g,%.12g,%.12g)-(%.12g,%.12g,%.12g)" % (e[0].collect { |c| c * @dbu } + e[1].collect { |c| c * @dbu })
  end

end

