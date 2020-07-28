CC = g++

build:  SimBox.cpp
		$(CC) -g -o run SimBox.cpp myrandom.cpp AABB.cc    
