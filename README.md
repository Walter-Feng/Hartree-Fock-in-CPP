# Hartree-Fock-in-CPP

=============================================================================

Three years ago, there was a Hartree Fock program written in C (https://github.com/Walter-Feng/Hartree-Fock). During that time, I only learned C, so I needed to write all the structures myself - tons of pointers, stupid recursion, and so on. The code was super inefficient, buggy, maybe with some memory leaks.

Things changed. I learned C++ and software engineering practice, read some good code, wrote a lot more code in different languages. Externally, C is no longer popular among undergraduates. Therefore, reviewing the original C-style Hartree-Fock code, I thought instantly "hey this code need a complete reconstruction", so making it a lot easier for newcomers to understand the math, run the code, and get a taste of what good code is.

"Well... but you could have used Python which is more friendly"

No. I refuse to use python for this project. Learn some static programming!

=============================================================================

Back to the topic, this is a Hartree-Fock program that aims to show how Hartree-Fock really works in every process, but in full C++ standard. (Hopefully) it is well documented with good variable names so that you can read the code without much effort and understand what is going on. A good tex document that illustrates the mathematics and C++ features is (probably) on the way.

As a comparison, the features in this program are:

   1. The programming style is completely in C++ - OOP, templates, but no more pointers.

   2. Tons of dependencies, including 
        Boost (you need to install yourself)
        
        Armadillo (same as above)
        
        GSL - GNU Scientific Library (same as above)
        
        sunqm/libcint (but only using its rys quadrature roots generator)
        
        catchorg/Catch2
        
        Taywee/args
        
        fmtlib/fmt
        
        nlohmann/json

   3. It is (in principle) able to use all the bases that have been published in Basis Set Exchange (https://www.basissetexchange.org).

   4. No special trick is used generally, but some C++ tricks may need to be learned.

   5. Templates may prevent understanding the code via some debug techniques.

And something good:

   1. GOOOOOOD PROGRAMMING STYLE

   2. The output is mostly handled by printer functions, so easier to mess around.

To clone the repository, don't forget to add "--recursive" to include all the dependencies,
```
git clone --recursive https://github.com/Walter-Feng/Hartree-Fock-in-CPP.git
```

CMake is used to generate the executable:
```
Hartree-Fock-in-CPP-repo-directory: $ mkdir build 
Hartree-Fock-in-CPP-repo-directory/build: $ cmake ..
Hartree-Fock-in-CPP-repo-directory/build: $ make
```

Type
```
./hfincpp -h
```
to unlock all options.

To run the example:
```
Hartree-Fock-in-CPP-repo-directory/build: $ ./hfincpp ../example/H2_6_31g.json
```

The output should be something like:
```
============================================================
|Iter|        Time / s |   Energy / a.u. |     Energy Diff |
============================================================
     0                 0       -0.90365255                 0
     1        2.0781e-05       -0.91596619      -0.012313641
     2        3.5345e-05       -0.91626362    -0.00029743166
     3        4.6462e-05       -0.91627108    -7.4536724e-06
     4        5.6896e-05       -0.91627126    -1.8785583e-07
============================================================

Total time elapsed: 0.032478466 s
```

Json file format is used for input, and due to lack of default parameters, any missing keywords will generate a non-trackable error. So for now it is advised to only change the values in example json input file.

Now the code also supports calculating gradient and perform a BFGS optimization. An example of optimizing H2O is shown below:

```
Performing SCF ...
============================================================
|Iter|        Time / s |   Energy / a.u. |     Energy Diff |
============================================================
     0                 0        -73.320014                 0
     1        3.8087e-05        -74.813641        -1.4936266
     2        6.4239e-05        -74.896802      -0.083161128
     3        8.6848e-05        -74.902545     -0.0057425592
     4       0.000108737        -74.903398    -0.00085350827
     5       0.000129691        -74.903606    -0.00020763838
     6       0.000149666         -74.90366    -5.3901557e-05
     7       0.000169352        -74.903674    -1.4139079e-05
     8       0.000188594        -74.903677    -3.7127233e-06
     9       0.000207793        -74.903678    -9.7514321e-07
============================================================

Calculating gradient ...
============================================================
|Atom| X / (Ha / bohr) | Y / (Ha / bohr) | Z / (Ha / bohr) |
============================================================
     H        -0.1256186     0.00088870278                 0
     O    -1.5781543e-12     -0.0017774056                 0
     H         0.1256186     0.00088870278                 0
============================================================

Performing optmization ...
============================================================
|Iter|        Time / s |   Energy / a.u. |     Grad / a.u. |
============================================================
     0                 0        -74.903678         0.1256186
     1       0.084930816        -74.920302        0.10866145
     2        0.17214759        -74.944935        0.06184493
     3        0.72308432        -74.953628       0.093235889
     4         0.8116638        -74.965104       0.015954198
     5         1.2036186        -74.965662        0.01045005
     6         1.2902973        -74.965866      0.0038400685
     7         1.7879789          -74.9659     0.00072399627
============================================================

Optimized geometry:
============================================================
|Atom|    X / Angstrom |    Y / Angstrom |    Z / Angstrom |
============================================================
     H        0.45780056      -0.086193387                 0
     O          1.889726         1.1172498                 0
     H         3.3216514      -0.086193387                 0
============================================================
Total time elapsed: 1.989305952 s

Process finished with exit code 0
```



   
