########reading nodes and elements########

gfx read nodes GENERICEXAMPLE NODE_DATA_MOONEY_RIVLIN.part0.exnode

gfx read elements GENERICEXAMPLE ELEMENT_DATA_MOONEY_RIVLIN.part0.exelem

########defining faces########

gfx define faces

########defining depedent fields##########

gfx define field Deformed_xx component DEPENDENT_FIELD.1  ;

gfx define field Deformed_yy component DEPENDENT_FIELD.2  ;

gfx define field Deformed_zz component DEPENDENT_FIELD.3  ;

gfx define field disp_xx add fields DEPENDENT_FIELD.1 GEOMETRY.x scale_factors 1 -1;

gfx define field disp_yy add fields DEPENDENT_FIELD.2 GEOMETRY.y scale_factors 1 -1;

gfx define field disp_zz add fields DEPENDENT_FIELD.3 GEOMETRY.z scale_factors 1 -1;

gfx define field DeformedGeometry component DEPENDENT_FIELD.1 DEPENDENT_FIELD.2 DEPENDENT_FIELD.3;


######## opening window to visualize results and introducing coordinate axis##########

gfx cre win 1 ;


gfx modify g_element "/" point glyph axes_xyz font BIG general size "2*2*2" select_on material black selected_material black ;


########creating creating a spectrum to display DEPENDENT_FIELD fields ##########

gfx create spectrum displacement_spectrum


######## Creating surface ########

gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material black data disp_xx spectrum displacement_spectrum selected_material default_selected render_shaded;

gfx modify spectrum displacement_spectrum autorange

########creating nodes with labels########

gfx modify g_element "/" node_points subgroup REGION coordinate DeformedGeometry LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label DeformedGeometry select_on material black selected_material default_selected;

######## Creating another surface to show undeformed GEOMETRY ########

gfx modify g_element "/" surfaces coordinate GEOMETRY tessellation default LOCAL select_on material black  spectrum displacement_spectrum selected_material default_selected render_wireframe;

########displaying color bar on the visualization window. ##########

gfx create colour_bar spectrum displacement_spectrum label_material black

gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on normalised_window_fit_left;


gfx modify window 1 image view_all

gfx modify window 1 background colour 1 1 1

