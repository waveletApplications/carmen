# Makefile

#-----------------------------------------------------
# CHOOSE THE COMPILER BY COMMENTING ONE OF THESE LINES
#-----------------------------------------------------

# Compiler for the sequential version
#CPP=g++

# Compiler for the parallel version
CPP=mpiCC

#-----------------------------------------------------
# COMPILER FLAGS
#-----------------------------------------------------

# Optimization flags
OFLAG=-O2

# Debug flags
DFLAG=-g3

#-----------------------------------------------------
# COMMENT STATIC FLAG FOR PARALLEL COMPILATOR mpiCC
#-----------------------------------------------------

# Static flag
SFLAG=-static

# Warning flags
WFLAG=-Wall -Wno-deprecated

#FLAGS=$(OFLAG) $(SFLAG) $(WFLAG)
#FLAGS=$(DFLAG) $(SFLAG) $(WFLAG)
FLAGS=$(WFLAG)

OBJECTS=\
AdaptTimeStep.o \
AddErrorL1.o \
AddErrorL2.o \
Advection.o \
AnalyticAverage.o \
Backup.o \
BC.o \
BoundaryRegion.o \
Burgers.o \
Cell.o \
ComputedTolerance.o \
CPUTimeRef.o \
CurlVelocity.o \
DigitNumber.o \
FileWrite.o \
FineMesh.o \
Flux.o \
GetBoundaryCells.o \
InitAverage.o \
InitTimeStep.o \
InitTree.o \
Limiter.o \
Matrix.o \
MinAbs.o \
MoveGauss.o \
Node.o \
NormMaxQuantities.o \
Parallel.o \
Parameters.o \
Performance.o \
PrintGrid.o \
PrintIntegral.o \
RadiationRate.o \
ReactionRate.o \
RefreshTree.o \
Remesh.o \
SchemeAUSM.o \
SchemeCentered.o \
SchemeENO2.o \
SchemeENO3.o \
SchemeMcCormack.o \
SchemeOsmp3Mc.o \
SchemeOsmp5Mc.o \
SchemeRoe.o \
ShowTime.o \
Sign.o \
SignJacobian.o \
Smagorinsky.o \
Source.o \
Step.o \
Sutherland.o \
TimeAverageGrid.o \
TimeEvolution.o \
Timer.o \
Vector.o \
View.o \
ViewEvery.o \
ViewIteration.o \
ViscousFlux.o \
main.o

HEADERS=\
Carmen.h \
FineMesh.h \
Node.h \
PreProcessor.h \
TimeAverageGrid.h \
Vector.h \
Cell.h \
Matrix.h \
Parameters.h \
PrintGrid.h \
Timer.h

# Compile source files
.cpp.o:
	$(CPP) -c $(FLAGS) $<

# Link objects
carmen: $(OBJECTS) $(HEADERS)
	$(CPP) -o carmen $(FLAGS) $(OBJECTS)

clean:
	rm -f *.o
	rm -f *~
