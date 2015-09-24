
THE_SOURCE = ismrmrd_manuscript_mrm.tex
RESPONSE = ismrmrd_response
THE_GOAL = ismrmrd_paper

SHA1 = $(shell git rev-parse HEAD)
$(shell echo $(SHA1)>githash.txt)

# Run twice to get all references right
ALL: $(RESPONSE) $(THE_GOAL)

$(RESPONSE):
	pdflatex  -jobname $(RESPONSE) "\def\isresponse{1} \input{$(THE_SOURCE)}"
	bibtex $(RESPONSE)
	pdflatex  -jobname $(RESPONSE) "\def\isresponse{1} \input{$(THE_SOURCE)}"
	pdflatex  -jobname $(RESPONSE) "\def\isresponse{1} \input{$(THE_SOURCE)}"

$(THE_GOAL):
	pdflatex  -jobname $(THE_GOAL) $(THE_SOURCE)
	bibtex $(THE_GOAL)
	pdflatex  -jobname $(THE_GOAL) $(THE_SOURCE)
	pdflatex  -jobname $(THE_GOAL) $(THE_SOURCE)

clean:
	rm -f *~ *.bbl *.blg *.pdf *.log *.aux *.bak *.dvi *.sav *.out \
	      *.nav *.snm *.toc *.fff *.lof *.tdo *.lox *.brf *.soc *.gz \
          githash.txt

