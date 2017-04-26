# read in description
gfx read node Cantilever.part0.exnode
gfx read element Cantilever.part0.exelem

# define deformed geometry
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3

gfx create window 1

# display deformed geometry
gfx define faces egroup "Region"
gfx modify g_element "Region" lines coordinate deformed_geom select_on material green selected_material default_selected
gfx modify g_element "Region" node_points coordinate deformed_geom glyph sphere General size "2*2*2" centre 0,0,0 font default select_on material gold selected_material default_selected

# display undeformed lines
gfx modify g_element "Region 1" lines select_on material green selected_material default_selected

gfx create axes length 5 material default
gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2
