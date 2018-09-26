# Project: SuperellipsoidS
# Makefile created by Dev-C++ 4.9.9.2

CPP  = g++.exe
CC   = gcc.exe
WINDRES = windres.exe
RES  = 
OBJ  = frustum.o funcoes.o arquivos.o lighting.o algebra.o draw.o geometry.o main.o text.o $(RES)
LINKOBJ  = frustum.o funcoes.o arquivos.o lighting.o algebra.o draw.o geometry.o main.o text.o $(RES)
LIBS =  -L"C:/Dev-Cpp/lib" -lglut32 -lglu32 -lopengl32 -lwinmm -lgdi32  -s -m3dnow 
INCS =  -I"C:/Dev-Cpp/include" 
CXXINCS =  -I"C:/Dev-Cpp/lib/gcc/mingw32/3.4.2/include"  -I"C:/Dev-Cpp/include/c++/3.4.2/backward"  -I"C:/Dev-Cpp/include/c++/3.4.2/mingw32"  -I"C:/Dev-Cpp/include/c++/3.4.2"  -I"C:/Dev-Cpp/include" 
BIN  = Superellipsoids.exe
CXXFLAGS = $(CXXINCS)   -m3dnow
CFLAGS = $(INCS)   -m3dnow
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before Superellipsoids.exe all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "Superellipsoids.exe" $(LIBS)

frustum.o: frustum.cpp
	$(CPP) -c frustum.cpp -o frustum.o $(CXXFLAGS)

funcoes.o: funcoes.cpp
	$(CPP) -c funcoes.cpp -o funcoes.o $(CXXFLAGS)

arquivos.o: arquivos.cpp
	$(CPP) -c arquivos.cpp -o arquivos.o $(CXXFLAGS)

lighting.o: lighting.cpp
	$(CPP) -c lighting.cpp -o lighting.o $(CXXFLAGS)

algebra.o: algebra.c
	$(CPP) -c algebra.c -o algebra.o $(CXXFLAGS)

draw.o: draw.c
	$(CPP) -c draw.c -o draw.o $(CXXFLAGS)

geometry.o: geometry.cpp
	$(CPP) -c geometry.cpp -o geometry.o $(CXXFLAGS)

main.o: main.cpp
	$(CPP) -c main.cpp -o main.o $(CXXFLAGS)

text.o: text.c
	$(CPP) -c text.c -o text.o $(CXXFLAGS)
