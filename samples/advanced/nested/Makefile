AMGDIR=../../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_linsolve -lamg_prec -lpsb_prec
EXTRA_LIBS=-lblas
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

STOKES_OBJS=amg_d_nested_block_prec_mod.o amg_d_nest_stokes_file_test.o
OPTIMAL_OBJS=amg_d_nested_block_prec_mod.o amg_d_nest_optimal_control_test.o

EXEDIR=./runs

all: runsd amg_d_nest_stokes_file_test amg_d_nest_optimal_control_test

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

amg_d_nest_stokes_file_test: $(STOKES_OBJS)
	$(FLINK) $(LINKOPT) $(STOKES_OBJS) -o amg_d_nest_stokes_file_test \
	$(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS) $(EXTRA_LIBS)
	/bin/mv amg_d_nest_stokes_file_test $(EXEDIR)

amg_d_nest_optimal_control_test: $(OPTIMAL_OBJS)
	$(FLINK) $(LINKOPT) $(OPTIMAL_OBJS) -o amg_d_nest_optimal_control_test \
	$(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS) $(EXTRA_LIBS)
	/bin/mv amg_d_nest_optimal_control_test $(EXEDIR)

clean:
	/bin/rm -f $(STOKES_OBJS) $(OPTIMAL_OBJS) *$(.mod) \
	$(EXEDIR)/amg_d_nest_stokes_file_test \
	$(EXEDIR)/amg_d_nest_optimal_control_test

lib:
	(cd ../../..; make library)

verycleanlib:
	(cd ../../..; make veryclean)
