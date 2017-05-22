#Read in the sequence of nodal positions.
for $i (1..201)
  {
	 $filename = sprintf("ActiveStrain_TransIso_%d.part0.exnode", $i);
	 
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i-1;
  }

#Read in the element description
gfx read elements ActiveStrain_TransIso.part0.exelem;

# define deformed geometry
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3
gfx define field "Fx_mag" component "del U_del n.1"
gfx define field "Fy_mag" component "del U_del n.2"
gfx define field "Fz_mag" component "del U_del n.3" 
gfx define field "Fa_vec" component "del U_del n.1" "del U_del n.2" "del U_del n.3"
gfx define field "Fa_mag" magnitude field "del U_del n"
gfx define field "Muscle" component Material.10

gfx create window 1

# display deformed geometry
gfx define faces egroup "Region"
gfx modify g_element "Region" lines coordinate deformed_geom select_on material default selected_material default_selected
gfx modify g_element "Region" node_points coordinate deformed_geom glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected


#gfx create axes length 5 material default
#gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2

gfx create time_editor

gfx modify g_element "Region" general clear circle_discretization 6 default_coordinate Geometry element_discretization "4*4*4" native_discretization none;
gfx modify g_element "Region" lines select_on material default data Muscle selected_material default_selected;
gfx modify g_element "Region" lines coordinate deformed_geom select_on material gold selected_material default_selected;
gfx modify g_element "Region" node_points coordinate deformed_geom glyph arrow_line general size "0.0*0.0*0.0" centre 0,0,0 font default variable_scale Fa_mag scale_factors "1*1*1" orientation Fa_vec select_on material gold selected_material default_selected;
gfx modify g_element "Region" element_points coordinate deformed_geom glyph arrow_line general size "10*1*1" orientation Fibre material green variable_scale Muscle discretization 5
gfx modify g_element "Region" element_points coordinate deformed_geom glyph arrow_line general size "10*1*1" orientation Fa_vec material red variable_scale Fa_mag discretization 5

gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 721 height 815;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 437.156 -414.189 -1684.26 interest_point 120.45 164.215 -1510.98 up_vector 0.180983 -0.189986 0.964961 view_angle 40 near_clipping_plane 6.81821 far_clipping_plane 2436.59 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines antialias 2 depth_of_field 0.0 fast_transparency blend_normal;
