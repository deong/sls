CXX=g++
CXXFLAGS=-Wall -O3 -funroll-loops -pipe
#CXXFLAGS=-Wall -ggdb -pipe

OBJ=configuration.o\
	factory.o\
	keywords.o\
	mtrandom.o\
	problems.o\
	strategy.o\
	utilities.o

SRC=chc.cpp\
	chromosome.cpp\
	comparator.cpp\
	convergence.cpp\
	crossover.cpp\
	ea.cpp\
	emoea.cpp\
	encoding.cpp\
	genitor.cpp\
	landscape.cpp\
	localsearch.cpp\
	metrics.cpp\
	morrls.cpp\
	move.cpp\
	mutation.cpp\
	neighborhood.cpp\
	nsga2.cpp\
	pfront.cpp\
	population.cpp\
	repair.cpp\
	replacement.cpp\
	rrls.cpp\
	sandm.cpp\
	selection.cpp\
	simanneal.cpp\
	simplega.cpp\
	sls.cpp\
	spea2.cpp\
	tabu.cpp\
	terminators.cpp\
	twophasels.cpp

sls : libsls.a ${SRC} slsmain.cpp
	${CXX} ${CXXFLAGS} -L. slsmain.cpp -o sls -lsls -lm

doc : sls.doxygen
	doxygen sls.doxygen

gapgen : gapgen.o libsls.a
	${CXX} ${CXXFLAGS} -L. gapgen.o -o gapgen -lsls -lm

ksgen : ksgen.o libsls.a
	${CXX} ${CXXFLAGS} -L. ksgen.o -o ksgen -lsls -lm

enumerate : libsls.a ${SRC} enumerate.cpp
	${CXX} ${CXXFLAGS} -L. enumerate.cpp -o enumerate -lsls -lm

mergepf : libsls.a ${SRC} mergepf.cpp
	${CXX} ${CXXFLAGS} -L. mergepf.cpp -o mergepf -lsls -lm

pearson : pearson.cpp utilities.o
	${CXX} ${CXXFLAGS} pearson.cpp utilities.o -o pearson -lm

test : libsls.a test.cpp comparator.cpp
	${CXX} ${CXXFLAGS} -L. test.cpp -o test -lsls -lm

libsls.a : ${OBJ}
	ar rcs libsls.a ${OBJ}

clean :
	rm -f *.o

distclean :
	rm -f *.o libsls.a sls mergepf enumerate pearson gapgen ksgen
