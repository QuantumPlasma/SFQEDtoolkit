In order to use SFQEDtoolkit into your Montecarlo or PIC code just follow the steps listed hereunder. SFQEDtoolkit can be used by both cpp and fortran codes, thus the instructions below contain between brackets the programming language which they refer to. 

1.  Copy the folders `coefficients`, `src_SFQEDtoolkit` (required for **cpp** and **fortran**) and `src_fortran_wrapper` (only for **fortran**), pasting them to the directory where your main code is.

2. Change your code using the features of the toolkit. The library is thought to be used only through the functions that are made available in the `SFQEDtoolkit_Interface.hpp` header file (for **cpp**) or in the `SFQEDtoolkit_Interface.f90` source file (for **fortran**). However you are completely free to bypass this "self-imposed restriction" and use the library as you please. Examples of its usage are shown inside the `example_cpp` (**cpp**) or `example_fortran` (**fortran**) folder provided with the library (so please go have a look!).

3. Compile your code and SFQEDtoolkit together: (**cpp**) -> if you have a `makefile` this is easily done by adjusting it (so to include the content of the `src_SFQEDtoolkit` folder to your compilation path) and running `make`, otherwise you can always proceed in the old-fashioned way, by compiling and linking everything through the command line; (**fortran**) -> same as in the **cpp** case. Additionally, it is important to compile the content of the `src_fortran_wrapper` directory before compiling the files where the library is used.

4. Give SFQEDtoolkit a try, and apply the instructions above to the examples (`example_cpp` and `example_fortran`) we provided.

