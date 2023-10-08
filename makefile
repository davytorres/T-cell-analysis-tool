FC=gfortran -O
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c

SOURCES = tcjnew.f readcellsnew.f observedtcell.f identify_tracks.f pearsonr.f median.f sort.f 

tcell: $(subst .f,.o,$(SOURCES))
	$(FC) -o $@ $+

clean:
	-rm -f *.o *.mod *.smod main
