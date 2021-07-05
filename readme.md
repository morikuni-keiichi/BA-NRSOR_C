# C codes for the BA-GMRES method preconditioned by NR-SOR inner iterations 

by Keiichi Morikuni and Ken Hayami

## License

This software is currently released under the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html).

If you use this code in research for publication, please cite the papers

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Issue 1, pages 1-22, 2013. DOI: [10.1137/110828472](https://doi.org/10.1137/110828472)
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Issue 1, pages 225-250, 2015. DOI: [10.1137/130946009](https://doi.org/10.1137/130946009)

For commercial use, please make a contact to
Keiichi Morikuni [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp).

## Usage

To compile the codes, execute the following:

```
$ make
```

Change the Fortran compiler given in the Makefile file if necessary.

To simply run the program with the default values of parameters on a test matrix RANDL7, execute the following:

```
$ ./main 
```

Then the program outputs the approximate solution data solution.dat, the result data info.dat, and the relative residual norm history data reshis.dat.  
A test matrix called RANDL7 in the compressed column storage (CCS) format is given in directory RANDL7. 
In addition, some spcific data is output in log.csv.
Parameters such as the number of NR-SOR inner iterations and the NR-SOR relaxation parameter are written in prm.dat.
The NR-SOR inner iteration parameters can be automatically tuned at each restart.

## Contacts

Please provide feedback to [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp) if you have any questions or suggestions.

Keiichi Morikuni, Ph.D.  

Affiliation: Faculty of Engineering, Information and Systems, University of Tsukuba  
Postal Address: 1-1-1 Tennodai, Tsukuba, Ibaraki 305-8573, Japan

Homepage URL: [http://researchmap.jp/KeiichiMorikuni/](http://researchmap.jp/KeiichiMorikuni/)

Supported by 
 The Graduate University for Advanced Studies (SOKENDAI)
 Shonan Village, Hayama, Kanagawa 240-0193 Japan

## Support

The Graduate University for Advanced Studies (SOKENDAI), Shonan Village, Hayama, Kanagawa 240-0193 Japan


### References

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Issue 1, pages 1-22, 2013. DOI: [10.1137/110828472](https://doi.org/10.1137/110828472)
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Issue 1, pages 225-250, 2015. DOI: [10.1137/130946009](https://doi.org/10.1137/130946009)
