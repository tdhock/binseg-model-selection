article.pdf: article.tex refs.bib
	pdflatex article
	bibtex article
	pdflatex article
	pdflatex article
article.tex: article.Rnw
	R -e "Sweave('article.Rnw')"
