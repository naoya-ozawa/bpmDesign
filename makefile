solidangle:	ssd_solidangle.cpp
	`root-config --cxx --cflags` -o solidangle ssd_solidangle.cpp `root-config --glibs`

trajectory:	ssd_trajectory.cpp
	`root-config --cxx --cflags` -o trajectory ssd_trajectory.cpp `root-config --glibs`


