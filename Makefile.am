SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4
if RUN_LINTER
lint: 
	$(LINTER) $(srcdir)/src/*.c
ALL_LINT_TARGET = lint
endif
