FC=gfortran -O
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c

SOURCES = tcjnew.f readcellsnew.f observedtcell.f identify_tracks.f calculate_volume_traversed.f calculate_mean_squared_displacement.f pearsonr.f median.f sort.f 

tcell: $(subst .f,.o,$(SOURCES))
	$(FC) -o $@ $+

clean:
	-rm -f *.o *.mod *.smod main
