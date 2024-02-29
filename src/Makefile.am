ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -Wall

lib_LTLIBRARIES = libpvcraze.la
libpvcraze_la_SOURCES=emalloc.c curve.c fitmodels.c Error.c quickselect.c curvemodel.c spline_iv.c splinefit.c KH_iv.c PC_iv.c diode_iv.c  bezier2_iv.c bezier3_iv.c pwl_iv.c apwl_iv.c neldermead.c diffevol.c furthestsampling.c lib-pv-craze.c tictoc.c emalloc.h curve.h fitmodels.h io.h Error.h quickselect.h curvemodel.h spline_iv.h splinefit.h KH_iv.h PC_iv.h diode_iv.h  bezier2_iv.h bezier3_iv.h pwl_iv.h apwl_iv.h neldermead.h diffevol.h furthestsampling.h lib-pv-craze.h  lib-pv-craze-structs.h tictoc.h  gen_ivanalyzeflags.sh register_models.sh collect_headers.sh libpvcraze.map
libpvcraze_la_LDFLAGS = -lm -Wl,--version-script=${srcdir}/libpvcraze.map -no-undefined  -version-number $(LIBVERSION)
include_HEADERS =  pvcraze.h

BUILT_SOURCES =  models.h model_table.h ivanalyzemessages.h ivanalyzeflags.h pvcraze.h Makefile_octave pyproject.toml setup.py
MOSTLYCLEANFILES = models.h model_table.h ivanalyzemessages.h ivanalyzeflags.h pvcraze.h octave-pv-craze.tar.gz setup.py python-pv-craze.tar.gz Makefile_octave pyproject.toml setup.py
modelfiles = bezier2_iv.h  KH_iv.h spline_iv.h bezier3_iv.h diode_iv.h PC_iv.h pwl_iv.h apwl_iv.h 

lib-pv-craze.$(OBJEXT): $(pv_craze_SOURCES) models.h model_table.h ivanalyzeflags.h
fitmodels.$(OBJEXT): models.h model_table.h fitmodels.c
curve.$(OBJEXT): ivanalyzemessages.h ivanalyzeflags.h curve.c curve.h
spline_iv.$(OBJEXT): splinefit.h

EXTRA_DIST = octave-pv-craze.cc Makefile_octave.in octave_configure.ac octave_makefile.am setup.py.in python-pv-craze.c pyproject.toml.in
dist_pkgdata_DATA = octave-pv-craze.tar.gz python-pv-craze.tar.gz

models.h model_table.h: $(modelfiles)
	${SHELL} ${srcdir}/register_models.sh  ${srcdir} $(modelfiles)
ivanalyzemessages.h ivanalyzeflags.h: curve.c gen_ivanalyzeflags.sh
	${SHELL} ${srcdir}/gen_ivanalyzeflags.sh  ${srcdir} curve.c
pvcraze.h: models.h
	${SHELL} ${srcdir}/collect_headers.sh ${srcdir}
Makefile_octave: Makefile_octave.in ivanalyzemessages.h ivanalyzeflags.h models.h model_table.h
	sed 's|version|$(subst :,.,$(VERSION))|g' ${srcdir}/Makefile_octave.in > Makefile_octave
	sed -i 's|cflags|$(CFLAGS)|g' Makefile_octave
	sed -i 's|lflags|$(LIBS)|g' Makefile_octave
	sed -i 's|srcdir|${srcdir}|g' Makefile_octave
pyproject.toml: pyproject.toml.in
	sed 's|__version__|$(subst :,.,$(VERSION))|g' ${srcdir}/pyproject.toml.in > pyproject.toml
setup.py: setup.py.in python-pv-craze.c
	sed 's|__version__|$(subst :,.,$(VERSION))|g' ${srcdir}/setup.py.in > setup.py
python-pv-craze.tar.gz: pyproject.toml setup.py
	mkdir -p python-pv-craze-$(subst :,.,$(VERSION))
	cp setup.py pyproject.toml ${srcdir}/python-pv-craze.c ${srcdir}/curve.h ${srcdir}/lib-pv-craze.h ${srcdir}/lib-pv-craze-structs.h pyproject.toml python-pv-craze-$(subst :,.,$(VERSION))
	tar -zcvf python-pv-craze.tar.gz python-pv-craze-$(subst :,.,$(VERSION))
	rm -rf  python-pv-craze-$(subst :,.,$(VERSION))
octave-pv-craze.tar.gz: Makefile_octave
	make -f Makefile_octave
	
