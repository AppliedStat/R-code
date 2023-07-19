rm fig.eps fig.ps fig-crop.p* 
latex  fig.tex
dvips -Ppdf -G0 fig.dvi 
ps2pdf fig.ps
pdfcrop fig.pdf
pdftops fig-crop.pdf
ps2eps -f fig-crop.ps 
