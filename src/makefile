CXX=g++
CXXFLAGS=-Wall -I. -Ibem5 -O3 -funroll-loops -pipe
#CXXFLAGS=-Wall -I. -Ibem5 -ggdb -pipe

OBJ=bem5/BEMLoop.o\
	bem5/Burfell.o\
	bem5/RotorProfiles.o\
	bem5/turbine.o\
	bem5/TurbineCost.o\
	bem5/WeibullAEP.o\
	factory.o\
	keywords.o\
	mtrandom.o\
	problems.o\
	strategy.o\
	utilities.o

SRC=bem5/BEMLoop.cpp\
	bem5/Burfell.cpp\
	bem5/RotorProfiles.cpp\
	bem5/turbine.cpp\
	bem5/TurbineCost.cpp\
	bem5/WeibullAEP.cpp\
	chc.cpp\
	chromosome.cpp\
	comparator.cpp\
	convergence.cpp\
	crossover.cpp\
	de.cpp\
	ea.cpp\
	emoea.cpp\
	encoding.cpp\
	genitor.cpp\
	landscape.cpp\
	localsearch.cpp\
	metrics.cpp\
	morrls.cpp\
	lsmove.cpp\
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

sls : libsls.a libkvp.a ${SRC} slsmain.cpp
	${CXX} ${CXXFLAGS} -L. slsmain.cpp -o sls -lsls -lkvp -lboost_regex -lm

doc : sls.doxygen
	doxygen sls.doxygen

gapgen : gapgen.o libsls.a
	${CXX} ${CXXFLAGS} -L. gapgen.o -o gapgen -lsls -lm

ksgen : ksgen.o libsls.a
	${CXX} ${CXXFLAGS} -L. ksgen.o -o ksgen -lsls -lm

enumerate : libsls.a libkvp.a ${SRC} enumerate.cpp
	${CXX} ${CXXFLAGS} -L. enumerate.cpp -o enumerate -lsls -lkvp -lboost_regex -lm

mergepf : libsls.a libkvp.a ${SRC} mergepf.cpp
	${CXX} ${CXXFLAGS} -L. mergepf.cpp -o mergepf -lsls -lkvp -lboost_regex -lm

pearson : pearson.cpp utilities.o
	${CXX} ${CXXFLAGS} pearson.cpp utilities.o -o pearson -lm

test : libsls.a test.cpp comparator.cpp
	${CXX} ${CXXFLAGS} -L. test.cpp -o test -lsls -lm

libsls.a : ${OBJ}
	ar rcs libsls.a ${OBJ}

libkvp.a :
	make -C kvparse libkvp.a
	cp kvparse/libkvp.a .

clean :
	rm -f *.o
	make -C kvparse clean

distclean :
	rm -f *.o libsls.a libkvp.a sls mergepf enumerate pearson gapgen ksgen
	make -C kvparse distclean
