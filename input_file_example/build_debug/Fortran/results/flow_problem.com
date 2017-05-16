########reading nodes and elements########

gfx read nodes GENERICEXAMPLE TimeStep_0.part0.exnode time 0


gfx read elements GENERICEXAMPLE ELEMENTDATA_FLUID_FLOW.part0.exelem

#### reading nodal values on different time steps ##################
for ($i = 1 ; $i < 11 ; $i++)
  {
	 gfx read nodes GENERICEXAMPLE TimeStep_$i.part0.exnode time $i;
  }

########defining faces########

gfx define faces

########defining depedent fields##########

gfx define field vel_x component DEPENDENTFIELD.1  ;

gfx define field vel_y component DEPENDENTFIELD.2  ;

gfx define field pressure component  DEPENDENTFIELD.3  ;

######## opening window to visualize results , witch to 2D view and introducing coordinate axis##########

gfx cre win 1 ;

gfx modify window 1 view parallel  eye_point -0.2 0.1 4 up_vector 0 1 0

gfx modify g_element "/" point glyph axes_xyz font BIG general size "2*2*2" select_on material black selected_material black ;


########creating surfaces with dependent fields ##########

gfx mod g_e REGION surfaces data vel_x spectrum default

########creating nodes with labels##########
#point label cmiss_number

gfx modify g_element REGION node_points glyph point  select_on material default selected_material default_selected;

########creating creating a spectrum to display dependent fields ##########

gfx create spectrum vel_spectrum autorange

gfx mod g_e REGION surfaces data vel_x spectrum vel_spectrum

gfx modify spectrum vel_spectrum autorange

########displaying color bar on the visualization window. ##########

gfx create colour_bar spectrum vel_spectrum label_material black

gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material black selected_material black normalised_window_fit_left;

################# creating time editor #####################################
gfx create time_editor

gfx timekeeper default stop set_time 0

gfx timekeeper default play speed 3



###########################################################################
gfx modify window 1 image view_all

gfx modify window 1 background colour 1 1 1

