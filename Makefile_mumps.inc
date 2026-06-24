# Begin orderings
LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

ISCOTCHDIR = $(topdir)/../../include
LSCOTCHDIR = $(topdir)/../../lib
ISCOTCH   = -I$(ISCOTCHDIR)
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

IMETISDIR = $(topdir)/../../include
LMETISDIR = $(topdir)/../../lib
IMETIS    = -I$(IMETISDIR)
LMETIS    = -L$(LMETISDIR) -lmetis

# Corresponding variables reused later
ORDERINGSF = -Dmetis -Dpord -Dscotch 
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(IMETIS) $(IPORD) $(ISCOTCH)
IORDERINGSC = $(IORDERINGSF)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
# -fPIC produces adress position independent code
CC = gcc -fPIC
FC = gfortran -fPIC
FL = gfortran -fPIC
AR = ar vr 
RANLIB = ranlib

LAPACK = -L$(topdir)/../lib/ -lopenblas
LIBBLAS = -L$(topdir)/../lib/ -lopenblas
LIBOTHERS = -lpthread

INCSEQ = -I$(topdir)/libseq
# libmpiseq linking for sequential build
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq/ -lmpiseq

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
# -fdefault-integer-8 force gfortran to use 8-Byte Integer (consistant with OpenBLAS, Scotch, METIS)
# -DINTSIZE64 stellt das C-Schnittstellen-Layout von MUMPS auf 8-Byte Integer um
# may add following FFlags to fixes https://github.com/coin-or-tools/ThirdParty-Mumps/issues/4 (only required for gcc/gfortran >= 10):
# -fallow-argument-mismatch -fallow-invalid-boz
OPTF    = -O3 -fopenmp -DBLR_MT -fallow-argument-mismatch -fdefault-integer-8
OPTL    = -O3 -fopenmp -DBLR_MT
OPTC    = -O3 -fopenmp -DBLR_MT -DINTSIZE64
#End Optimized options

INCS = $(INCSEQ)
LIBS = -L$(libdir) -lmumps_common $(LPORD) $(LSCOTCH) $(LMETIS) $(LIBSEQ) $(LAPACK) $(LIBBLAS)
LIBSEQNEEDED = libseqneeded
