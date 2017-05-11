FUKUFU
=======
Fukufu small program written in c++ to calculate the Fukui function and other reactivity indexes from Gaussian outputs.

USAGE
------
A Gaussian (preferably 09) single-point open-shell calculation is needed, run with "pop=full iop(3/33=4)" options.
Then run with:

```
./fukufu.o gaussian_log_file.log
```

An output containing several reactivity indexes will be created.


INSTALATION
------------
Just compile the code with the following command line:

```
g++ src.cpp -o fukufu.o -std=c++11
```

g++ may be replaced with you favourite c++ compiler.
