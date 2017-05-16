########reading nodes and elements########

gfx read nodes GENERICEXAMPLE NODE_DATA_LAPLACE.part0.exnode

gfx read elements GENERICEXAMPLE ELEMENT_DATA_LAPLACE.part0.exelem

########defining faces########

gfx define faces

########defining depedent fields##########

gfx define field Temperature component DEPENDENT.1  ;


######## opening window to visualize results , witch to 2D view and introducing coordinate axis##########

gfx cre win 1 ;


gfx modify g_element "/" node_points glyph axes_xyz font BIG general size "2*2*2" select_on material black selected_material black ;


########creating nodes with labels##########
gfx modify g_element "/" general clear;
gfx modify g_element "/" node_points subgroup REGION coordinate GEOMETRY LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label cmiss_number select_on material black selected_material default_selected;


########creating creating a spectrum to display temperature field ##########

gfx create spectrum temperature_spectrum autorange

########creating surfaces with Temperature field ##########

gfx modify g_element "/" surfaces coordinate GEOMETRY tessellation default LOCAL select_on material black data Temperature spectrum temperature_spectrum selected_material default_selected render_shaded;


gfx modify spectrum temperature_spectrum


########displaying color bar on the visualization window. ##########

gfx create colour_bar spectrum temperature_spectrum label_material black

gfx modify g_element "/" node_points glyph colour_bar general size "1*1*1" centre 0,0,0 select_on  normalised_window_fit_left;



gfx modify window 1 image view_all

gfx modify window 1 background colour 1 1 1

