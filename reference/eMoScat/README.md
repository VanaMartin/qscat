Quantum Scattering Computational (QSCAT) package
================================================

This document contains basic instruction for compiling the eMoScat
project.

PRE-REQUISITIES
---------------

standard packages:
openblas
        NOTE:
        compile and than run: $make PREFIX=/opt/OpenBlas install

third party commercial packages:
intel compiler

TODO list
---------
-   Redo all possible save/read binary methods via BinaryStreamInterface (DONE)
+   Build shallow copy semantics via reference counting
+   Clean code from TODOs
+   Finish the documentation of libs
+   Allow static and dynamic libraries compilation
+   Redesign the Input library
+   Allow export to json via python scripts
+   Redesign python wrapper (without boost)
+   FFT (or FBT, FCT) via matrix vector multiplication (new class Fourier Filter)
+   Function initialization as special costructor to gVector resp. derived class

References and Acknowledgements:
--------------------------------

+   coulcc.f
+   picojson.h
+   openblas


BUGS:
=====

+   Segfault on inconsistent LCP load!
+   Evolution FAIL!
