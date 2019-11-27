# Makefile to run pieces of the analysis, generate images
# and produce the final document.

PAPER=report

images:
	echo "MAKING IMAGES"

paper:
	pdflatex $(PAPER).tex
	bibtex $(PAPER).aux
	pdflatex $(PAPER).tex
	pdflatex $(PAPER).tex

clean:
	rm -f $(PAPER).aux $(PAPER).log \
        $(PAPER).pdf $(PAPER).blg $(PAPER).bbl ${PAPER}.spl \
        $(PAPER)*backup
