bin_PROGRAMS = quick_example position_table solid_angle compare_horizons
position_table_SOURCES = position_table.c
quick_example_SOURCES = quick_example.c
compare_horizons_SOURCES = compare_horizons.c
solid_angle_SOURCES = solid_angle.c
position_table_LDADD   = ../src/libaephem.la -lm
quick_example_LDADD   = ../src/libaephem.la -lm
solid_angle_LDADD   = ../src/libaephem.la -lm
compare_horizons_LDADD   = ../src/libaephem.la -lm
position_table_LDFLAGS = -Wl,--no-as-needed
quick_example_LDFLAGS = -Wl,--no-as-needed
solid_angle_LDFLAGS = -Wl,--no-as-needed
compare_horizons_LDFLAGS = -Wl,--no-as-needed
INCLUDES	 = -I@top_srcdir@/src
