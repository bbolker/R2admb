RM = rm -f
CP = cp
TOUCH = touch
REXE = R
ROPTS = --vanilla
RCMD = $(REXE) $(ROPTS) CMD
RSCRIPT = $(R_HOME)/bin/Rscript $(ROPTS)

CC = gcc
CXX = g++
CFLAGS = -O4 -Wall
CXXFLAGS = -O4 -Wall
LDLIBS =
FC = gfortran
FFLAGS = -O4
INSTALL = install

LATEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEIDX = makeindex
F2D = fig2dev

DVIPS = dvips -tletter -Ppdf
PSTOPDF = ps2pdf
EPSTOPDF = epstopdf
XMLCMD = mk4ht mzlatex
HTMLCMD = mk4ht htlatex
TTH = ./tthconv
TTM = ttmconv

default:

.PHONY:

%.tex: %.Rnw
	$(RSCRIPT) -e "library(weaver); Sweave(\"$*\",driver=weaver());"

%.R: %.Rnw
	$(RCMD) Stangle $*

%.html: %.tex
	$(HTMLCMD) $*

%.xml: %.tex
	$(XMLCMD) $*

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.ps
	$(PSTOPDF) $*.ps $@

%.ps: %.dvi
	$(DVIPS) $*.dvi -o $@

%.pdf: %.fig
	$(F2D) -Lpdf $*.fig $@

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

%.eps: %.fig
	$(F2D) -Leps $*.fig $@

%.png: %.fig
	$(F2D) -Lpng $*.fig $@

%.bbl: %.tex
	-$(LATEX) $*
	$(BIBTEX) $*

%.idx: %.tex
	-$(LATEX) $*

%.ind: %.idx
	$(MAKEIDX) $*

%: %.cc
	$(CXX) $(CXXFLAGS) $*.cc $(LDLIBS) -o $*

%: %.c
	$(CC) $(CFLAGS) $*.c $(LDLIBS) -o $*

%: %.f
	$(FC) $(FFLAGS) $*.f $(LDLIBS) -o $*

clean:
	$(RM) *.o *.tex *.log *.aux *.toc *.blg *.bbl *.out *-???.pdf Rplots.ps

.SECONDARY:
