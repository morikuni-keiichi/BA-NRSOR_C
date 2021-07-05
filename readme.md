# C codes for the BA-GMRES method preconditioned by the NR-SOR inner iterations 

by Keiichi Morikuni and Ken Hayami

## License

This software is currently released under the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html).

If you use this code in research for publication, please cite the papers

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Issue 1, pages 1-22, 2013. DOI: [10.1137/110828472](https://doi.org/10.1137/110828472)
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Issue 1, pages 225-250, 2015. DOI: [10.1137/130946009](https://doi.org/10.1137/130946009)

For commercial use, please make a contact to
Keiichi Morikuni [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp).

## Usage

Maintained by
 Keiichi Morikuni <morikuni@cs.cas.cz>
 Institute of Computer Science, Academy of Sciences of the Czech Republic
 Pod Vodarenskou vezi 271/2, 182 07 Prague 8, Czech Republic


Supported by 
 The Graduate University for Advanced Studies (SOKENDAI)
 Shonan Village, Hayama, Kanagawa 240-0193 Japan


This software is currently released under the GNU General Public License
http://www.gnu.org/copyleft/gpl.html


If you use this code in research for publication, please cite the papers

Morikuni, K. and Hayami, K., 
Inner-iteration Krylov subspace methods for least squares problems,
SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, 
pages 1-22, 2013.

Morikuni, K. and Hayami, K., 
Convergence of inner-iteration GMRES methods for least squares problems (Revised), 
NII Technical Report, National Institute of Informatics, NII-2013-004E 1-24, 2013.
http://www.nii.ac.jp/TechReports/13-004E.html


For the commercial use, please make a contact to 
Keiichi Morikuni morikuni@cs.cas.cz


BA-GMRES preconditioned by NR-SOR inner iterations involves the following files:

   main.c 		last update October 7, 2014
   solver.c		last update November 19, 2014
   sub.c      last update November 19, 2014
   func.c	
   globvar.c	last update October 7, 2014
   prm.dat		last update October 7, 2014
   plot.plt		(may be used if gnuplot is available)
   makefile		last update October 7, 2014
   header.h   last update October 7, 2014
   readme     last update November 19, 2014

A test matrix called RANDL7 in the compressed column storage (CCS) format 
is given in RANDL7 directory. Parameters such as the number of NR-SOR inner 
iterations and the NR-SOR relaxation parameter are written in prm.dat.
The NR-SOR inner iteration parameters can be automatically tuned at each restart.


To compile the codes and run the program, proceed as follows:

  $ make 
  $ ./main RANDL7/

Then the program outputs the approximate solution data solution.dat,
the result data info.dat, and the relative residual norm history data reshis.dat.  
In addition, some spcific data is output in log.csv.

Please provide feedback if you have any questions or suggestions.
morikuni@cs.cas.cz


References
[1] Morikuni, K. and Hayami, K., 
Inner-iteration Krylov subspace methods for least squares problems,
SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, pages 1-22, 2013.
[2] Morikuni, K. and Hayami, K.
Convergence of inner-iteration GMRES methods for least squares problems (Revised), 
NII Technical Report, National Institute of Informatics, NII-2013-004E 1-24, 2013.
http://www.nii.ac.jp/TechReports/13-004E.html