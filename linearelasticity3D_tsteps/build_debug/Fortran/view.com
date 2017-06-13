#Read in the sequence of nodal positions.
for $i (1..5)
  {
	 $filename = sprintf("UniaxialExtension2D_%d.part0.exnode", $i);
	 
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elements LinearElasticity2DPlaneStress.part0.exelem;

# define deformed geometry
gfx define field "deformed_geom" component Dependent.1 Dependent.2

gfx create window 1

# display deformed geometry
gfx define faces egroup "Region 1"
gfx modify g_element "Region 1" lines coordinate deformed_geom select_on material default selected_material default_selected
gfx modify g_element "Region 1" node_points coordinate deformed_geom glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

gfx create axes length 5 material default
gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2

gfx create time_editor

# ###
#gfx define field deformed_geom coordinate_system rectangular_cartesian composite Dependent;
#gfx modify spectrum default clear overwrite_colour;
#gfx modify spectrum default linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;
#gfx create material black normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 1 shininess 0.2;
#gfx create material blue normal_mode ambient 0 0 0.5 diffuse 0 0 1 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
#gfx create material bone normal_mode ambient 0.7 0.7 0.6 diffuse 0.9 0.9 0.7 emission 0 0 0 specular 0.1 0.1 0.1 alpha 1 shininess 0.2;
#gfx create material default normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
#gfx create material default_selected normal_mode ambient 1 0.2 0 diffuse 1 0.2 0 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
#gfx create material gold normal_mode ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
#gfx create material gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 1 shininess 0.2;
#gfx create material green normal_mode ambient 0 0.5 0 diffuse 0 1 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.1;
#gfx create material muscle normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 1 shininess 0.2;
#gfx create material red normal_mode ambient 0.5 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
#gfx create material silver normal_mode ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
#gfx create material tissue normal_mode ambient 0.9 0.7 0.5 diffuse 0.9 0.7 0.5 emission 0 0 0 specular 0.2 0.2 0.3 alpha 1 shininess 0.2;
#gfx create material transparent_gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 0 shininess 0.2;
#gfx create material white normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
#gfx modify g_element "Region 1" general clear circle_discretization 6 default_coordinate Geometry element_discretization "4*4*4" native_discretization none;
#gfx modify g_element "Region 1" lines select_on material default selected_material default_selected;
#gfx modify g_element "Region 1" lines coordinate deformed_geom select_on material gold selected_material default_selected;
#gfx modify g_element "Region 1" node_points coordinate deformed_geom glyph sphere general size "5*5*5" centre 0,0,0 font default select_on material gold selected_material default_selected;
#gfx modify g_element "Region 1" node_points coordinate deformed_geom glyph arrow_line general size "5*5*5" centre 0,0,0 font default label "del U_del n" orientation "del U_del n" scale_factors "0.01*0.01*0.01" select_on material gold data Dependent spectrum default selected_material default_selected;
#gfx create window 1 double_buffer;
#gfx modify window 1 image scene default light_model default;
#gfx modify window 1 image add_light default;
#gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25 width 586 height 570;
#gfx modify window 1 set current_pane 1;
#gfx modify window 1 background colour 0 0 0 texture none;
#gfx modify window 1 view parallel eye_point 60 80.3333 409.068 interest_point 60 80.3333 0 up_vector -0 1 -0 view_angle 40 near_clipping_plane 4.09068 far_clipping_plane 1461.87 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
#gfx modify window 1 overlay scene none;
#gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines antialias 2 depth_of_field 0.0 fast_transparency blend_normal;

