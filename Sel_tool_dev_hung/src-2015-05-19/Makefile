
f = genomic-selection
g = plabsim
h = SelectionTools

# Linux

MTLin: cpMTLin 
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h)-64.so
	rm  -f -f *.o 

MTLinGSLO: cpMTLinGSLO
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h)-64.so
	rm  -f -f *.o 

MTLinOB: cpMTLinOB
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h)-64.so
	rm  -f -f *.o 

MTLinMKL: cpMTLinMKL
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h)-64.so
	rm  -f -f *.o 


STLin: cpST
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h)-64.so
	rm  -f -f *.o 

# Windows

MTWin: cpMTWin 
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h).dll
	rm  -f -f *.o 


STWin: cpST
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h).dll
	rm  -f -f *.o 

STWinOB: cpSTWinOB
	rm  -f -f *.o 
	R CMD SHLIB $(f).c $(g).c -o $(h).dll
	rm  -f -f *.o 


# Copy Makevars

cpST:
	cp Makevars.ST    Makevars

cpSTWinOB:
	cp Makevars.STWinOB Makevars


cpMTLin:
	cp Makevars.MTLin Makevars

cpMTWin:
	cp Makevars.MTWin Makevars

cpMTLinGSLO:
	cp Makevars.MTLinGSLO Makevars

cpMTLinFlame:
	cp Makevars.MTLinFlame Makevars

cpMTLinOB:
	cp Makevars.MTLinOB Makevars

cpMTLinMKL:
	cp Makevars.MTLinMKL Makevars



.PHONY : clean
clean : 
	rm  -f -f *.o 
	rm  -f -f *~ 
	rm  -f -f #* 
