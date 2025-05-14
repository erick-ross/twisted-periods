# Twisted Periods

This repository contains code to supprt the paper [Twisted Periods of Modular Forms]().
In particular, we verified Conjectures 6.1 - 6.4 for `D,K <= ???`.
One can run the code by the following:

```
$ # Conjecture 6.1
$ sage compute.sage dets_l_a 50 50
$ # Conjecture 6.2
$ sage compute.sage dets_l_b 50 50
$ # Conjecture 6.3
$ sage compute.sage dets_chi_a 50 50
$ # Conjecture 6.4
$ sage compute.sage dets_chi_b 50 50
$ # Compute the a and b coefficients
$ sage compute.sage coeffs_a 16 15 > coeffs_a.txt
$ sage compute.sage coeffs_b 16 15 > coeffs_b.txt
```
