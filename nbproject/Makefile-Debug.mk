#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include TensorVisuAdvec-Sphere-Makefile.mk

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/sph.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/tensorHeap.o \
	${OBJECTDIR}/tensorglyph.o \
	${OBJECTDIR}/octree.o \
	${OBJECTDIR}/quicksort.o \
	${OBJECTDIR}/particle.o \
	${OBJECTDIR}/newGrid.o \
	${OBJECTDIR}/gcgTensor.o \
	${OBJECTDIR}/optimizedGrid.o \
	${OBJECTDIR}/squadric.o \
	${OBJECTDIR}/tensorgenerate.o \
	${OBJECTDIR}/superquadric.o \
	${OBJECTDIR}/tensorvisu.o \
	${OBJECTDIR}/tensorgrid.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-lglui -lglut -lGLU -lGL -ljpeg -lgcg64 -ljpeg

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/agoravai

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/agoravai: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/agoravai ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/sph.o: sph.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/sph.o sph.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/tensorHeap.o: tensorHeap.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/tensorHeap.o tensorHeap.cpp

${OBJECTDIR}/tensorglyph.o: tensorglyph.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/tensorglyph.o tensorglyph.cpp

${OBJECTDIR}/octree.o: octree.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/octree.o octree.cpp

${OBJECTDIR}/quicksort.o: quicksort.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/quicksort.o quicksort.cpp

${OBJECTDIR}/particle.o: particle.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/particle.o particle.cpp

${OBJECTDIR}/newGrid.o: newGrid.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/newGrid.o newGrid.cpp

${OBJECTDIR}/gcgTensor.o: gcgTensor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/gcgTensor.o gcgTensor.cpp

${OBJECTDIR}/optimizedGrid.o: optimizedGrid.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/optimizedGrid.o optimizedGrid.cpp

${OBJECTDIR}/squadric.o: squadric.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/squadric.o squadric.cpp

${OBJECTDIR}/tensorgenerate.o: tensorgenerate.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/tensorgenerate.o tensorgenerate.cpp

${OBJECTDIR}/superquadric.o: superquadric.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/superquadric.o superquadric.cpp

${OBJECTDIR}/tensorvisu.o: tensorvisu.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/tensorvisu.o tensorvisu.cpp

${OBJECTDIR}/tensorgrid.o: tensorgrid.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/GL -MMD -MP -MF $@.d -o ${OBJECTDIR}/tensorgrid.o tensorgrid.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/agoravai

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
