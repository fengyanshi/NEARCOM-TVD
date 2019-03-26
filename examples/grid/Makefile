#  Makefile for ripple, linux version
#  fyshi oct 07 2002

#  Build everything.

all:  ripple

#  Build ripple.

ripple: ripple.o implctp.o porousset.o apbp.o prtplt.o setup.o aset.o initvoff.o setnf.o bc.o vtilde.o tension.o vofadv.o function.o wavemk.o turbmut.o klmodel.o comdk1.h iccgdk.h
	gfortran -o ripple ripple.o implctp.o porousset.o apbp.o prtplt.o setup.o aset.o initvoff.o setnf.o bc.o vtilde.o tension.o vofadv.o function.o wavemk.o turbmut.o klmodel.o
ripple.o: ripple.f comdk1.h
	gfortran -c -O ripple.f
implctp.o: implctp.f comdk1.h iccgdk.h
	gfortran -c -O implctp.f
porousset.o: porousset.f comdk1.h
	gfortran -c -O porousset.f
apbp.o: apbp.f comdk1.h
	gfortran -c -O apbp.f
prtplt.o: prtplt.f comdk1.h
	gfortran -c -O prtplt.f
setup.o: setup.f comdk1.h
	gfortran -c -O setup.f
aset.o: aset.f comdk1.h
	gfortran -c -O aset.f
initvoff.o: initvoff.f comdk1.h
	gfortran -c -O initvoff.f
setnf.o: setnf.f comdk1.h
	gfortran -c -O setnf.f
bc.o: bc.f comdk1.h
	gfortran -c -O bc.f
vtilde.o: vtilde.f comdk1.h
	gfortran -c -O vtilde.f
tension.o: tension.f comdk1.h
	gfortran -c -O tension.f
vofadv.o: vofadv.f comdk1.h
	gfortran -c -O vofadv.f
function.o: function.f comdk1.h
	gfortran -c -O function.f
wavemk.o: wavemk.f comdk1.h
	gfortran -c -O wavemk.f
klmodel.o: klmodel.f comdk1.h
	gfortran -c -O klmodel.f
turbmut.o: turbmut.f comdk1.h
	gfortran -c -O turbmut.f

#  Clean up the directory.

clean:
	rm *.o
