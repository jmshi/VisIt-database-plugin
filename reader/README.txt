
	Reader to read in rst data dump from athena++ and export data in
	*.debug files for debug. Most of the code are striped from athena++
	directly.  
	
	Currently, no support for 
	(a) parallel I/O 
	(b) user defined mesh data. 
	
	Things haven't been tested yet: 
	(a) non-uniform grid 
	(b) rst files with magnetic field data

	To start,
	1) run  `make all'  to compile and get executable named `reader'
	2) run  `./reader -r rst_filename'    e.g., rst_filename = amr.00000.rst
	3) coordinates for each meshblocks are then dumped in x1f.debug,
	x2f.debug and x3f.debug; conservative variables are in con.debug,
	field in mag.debug
