#Read in the sequence of nodal positions.
for $i (1..50)
  {
	 $filename = sprintf("ActiveStrain_TransIso_%d.part0.exnode", $i);
	 
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
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
gfx modify g_element "Region" lines select_on material default selected_material default_selected;
gfx modify g_element "Region" lines coordinate deformed_geom select_on material gold selected_material default_selected;
gfx modify g_element "Region" node_points coordinate deformed_geom glyph arrow_line general size "0.0*0.0*0.0" centre 0,0,0 font default variable_scale Fa_mag scale_factors "10*10*10" orientation Fa_vec select_on material gold selected_material default_selected;
gfx modify g_element "Region" element_points coordinate deformed_geom glyph arrow_line general size 0.1 orientation Fibre material red	
