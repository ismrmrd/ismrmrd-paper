
CURRENT_VERSION = ismrmrd_manuscript_mrm
THE_GOAL	= ismrmrd_paper
TARGETS		= $(THE_GOAL)

$(THE_GOAL).pdf: $(CURRENT_VERSION).tex
	latex $(CURRENT_VERSION).tex
	dvipdf $(CURRENT_VERSION).dvi $(THE_GOAL).pdf

# so that when latex is re-run, it will get all references right

update:
	touch $(CURRENT_VERSION).tex

really:
	rm -f *.ps $(THE_GOAL).pdf $(CURRENT_VERSION).pdf

clean:
	rm -f *~ *.log *.aux *.bak *.dvi *.sav *.out *.nav *.snm *.toc *.fff *.lof

