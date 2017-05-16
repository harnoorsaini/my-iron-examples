#Read in the sequence of nodal positions.
for $i (1..50)
  {
	 $filename = sprintf("ActiveStrain_Isometric_%d.part0.exnode", $i);
	 
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elements ActiveStrain_Isometric.part0.exelem;


gfx define field Dependent coordinate_system rectangular_cartesian finite_element number_of_components 4 field real component_names 1 2 3 4;
gfx define field Fibre coordinate_system fibre finite_element number_of_components 3 anatomical real component_names 1 2 3;
gfx define field Geometry coordinate_system rectangular_cartesian finite_element number_of_components 3 coordinate real component_names x y z;
gfx define field Material coordinate_system rectangular_cartesian finite_element number_of_components 5 field real component_names 1 2 3 4 5;
# field Region created by other commands
# field Region.cmiss_mesh_1d created by other commands
# field Region.cmiss_mesh_2d created by other commands
# field Region.cmiss_mesh_3d created by other commands
# field Region.cmiss_nodes created by other commands
gfx define field cmiss_number coordinate_system rectangular_cartesian cmiss_number;
gfx define field deformed_geom coordinate_system rectangular_cartesian composite Dependent.1 Dependent.2 Dependent.3;
gfx define field "del U_del n" coordinate_system rectangular_cartesian finite_element number_of_components 4 field real component_names 1 2 3 4;
gfx define field xi coordinate_system rectangular_cartesian xi_coordinates;
gfx define field xxx coordinate_system rectangular_cartesian composite "del U_del n.1";
gfx define field xya coordinate_system rectangular_cartesian composite "del U_del n.1" "del U_del n.2" "del U_del n.3";
gfx define field Fa coordinate_system rectangular_cartesian composite "del U_del n.1" "del U_del n.2" "del U_del n.3";
gfx define field Fx coordinate_system rectangular_cartesian composite "del U_del n.1";
gfx define field Fy coordinate_system rectangular_cartesian composite "del U_del n.2";
gfx define field Fz coordinate_system rectangular_cartesian composite "del U_del n.3";
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx create material black normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 1 shininess 0.2;
gfx create material blue normal_mode ambient 0 0 0.5 diffuse 0 0 1 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material bone normal_mode ambient 0.7 0.7 0.6 diffuse 0.9 0.9 0.7 emission 0 0 0 specular 0.1 0.1 0.1 alpha 1 shininess 0.2;
gfx create material default normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material default_selected normal_mode ambient 1 0.2 0 diffuse 1 0.2 0 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material gold normal_mode ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 1 shininess 0.2;
gfx create material green normal_mode ambient 0 0.5 0 diffuse 0 1 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.1;
gfx create material muscle normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 1 shininess 0.2;
gfx create material red normal_mode ambient 0.5 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material silver normal_mode ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material tissue normal_mode ambient 0.9 0.7 0.5 diffuse 0.9 0.7 0.5 emission 0 0 0 specular 0.2 0.2 0.3 alpha 1 shininess 0.2;
gfx create material transparent_gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 0 shininess 0.2;
gfx create material white normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1023 height 780;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 4.28141 -0.921018 1.48454 interest_point 0.796377 0.463125 0.466735 up_vector -0.111444 0.391166 0.913548 view_angle 65.3828 near_clipping_plane 0.0388552 far_clipping_plane 13.8855 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines antialias 2 depth_of_field 0.0 fast_transparency blend_normal;