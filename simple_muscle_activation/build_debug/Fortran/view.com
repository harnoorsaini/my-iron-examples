#Read in the sequence of nodal positions.
for $i (1..50)
  {
	 $filename = sprintf("LargeUniaxialExtension_%d.part0.exnode", $i);
	 
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elements LargeUniaxialExtension.part0.exelem;

# define deformed geometry
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3

gfx create window 1

# display deformed geometry
gfx define faces egroup "Region"
gfx modify g_element "Region" lines coordinate deformed_geom select_on material default selected_material default_selected
gfx modify g_element "Region" node_points coordinate deformed_geom glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

gfx create axes length 5 material default
gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2

gfx create time_editor

