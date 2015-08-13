
CURRENT_VERSION = ismrmrd_manuscript_mrm
THE_GOAL	= ismrmrd_paper
TARGETS		= $(THE_GOAL)

SHA1 = $(shell git rev-parse HEAD)


$(THE_GOAL).pdf: $(CURRENT_VERSION)_hash.tex
	latex $(CURRENT_VERSION)_hash.tex
	bibtex $(CURRENT_VERSION)_hash
	latex $(CURRENT_VERSION)_hash.tex
	latex $(CURRENT_VERSION)_hash.tex
	dvipdf $(CURRENT_VERSION)_hash.dvi $(THE_GOAL).pdf

$(CURRENT_VERSION)_hash.tex: $(CURRENT_VERSION).tex
	sed -e s/MANUSCRIPT_SHA1/$(SHA1)/ $(CURRENT_VERSION).tex > $(CURRENT_VERSION)_hash.tex

# so that when latex is re-run, it will get all references right

update:
	touch $(CURRENT_VERSION)_hash.tex

really:
	rm -f *.ps $(THE_GOAL).pdf $(CURRENT_VERSION).pdf

clean:
	rm -f *~ *.bbl *.blg *.pdf *.log *.aux *.bak *.dvi *.sav *.out *.nav *.snm *.toc *.fff *.lof $(CURRENT_VERSION)_hash.tex

