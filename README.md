In order to include and compile SFQEDtoolkit into your Montecarlo or PIC code, you first have to export the path to SFQEDtoolkit's source directories, that is "src_SFQEDtoolkit_cpp" and "src_SFQEDtoolkit_fortran", into the environment variable SFQED_TOOLKIT_USER. Thus, we suggest you to open the .bash or .bash_profile file, in your home directory, and add this line

```export SFQED_TOOLKIT=location/of/SFQEDtoolkit```

and restart your session or simply source your .bash_profile with

```source .bash_profile```

Next, follow one of the following options, according to the language your main code is written with.

- cpp:
    in principle you should be able to successfully compile the SFQEDtoolkit library by simply copying and pasting the "src_SFQEDtoolkit_cpp" folder, which is bundled with this SFQEDtoolkit's repository you downloaded from github, into the directory that corresponds to your compilation path. In this way, you are perfectly allowed to launch your usual makefile. The library is thought to be used only through the functions that are made available in the "SFQEDtoolkit_Interface.hpp" header file, an example of its usage is shown inside the example_cpp folder provided with the library (so please go have a look!), but you are completely free to bypass this "self-imposed restriction" and use the library as you please. If you choose to ignore the user interface introduced with "SFQEDtoolkit_Interface.hpp", then bear in mind that the main methods of the library are defined in the class SFQED_Processes ("SFQED_Processes.h").
    In case you do not have any makefile you can always use the one available in the "example_cpp" directory. If you choose this option: open the example_cpp/makefile file and switch the variable "YOUR_COMPILATION_PATH" to the path where your code is stored, save and run 
    ```make```
    **But be careful**: it is very likely that you will need to compile your program with very specific compilator and/or linker flags (like those needed for openmp, for instance), so make sure to add those flags (it should not be too complicated). After the compilation the executable "cpp_tester" is produced.

- fortran:
    in this case things get a little tricky. First of all, the toolkit can be easily mounted on your fortran code by employing the "use SFQEDtoolkit_Interface" statement. This allows you to use all the functions defined in the "SFQEDtoolkit_Interface.f90" file, provided inside the folder "src_SFQEDtoolkit_fortran" (please check the "example_fortran" directory to see a demonstration of its usage).
    Regarding the compilation instead, we need to produce a hybrid f.90/cpp executable. This is achieved by first compiling all the source files (both .f90 and .cpp) into .o object file using the -c option of the correspoding compiler. Needless to say, during this phase, if the files you are compiling are in different folders, remember to specify explicitly the path to those through the include -I/paths option. After the .o files have been created you need to link them to produce the actual executable. Inside the "example_fortran" directory we provided a makefile: if you open that file and change the "YOUR_COMPILATION_PATH" variable to the path where your fortran code is stored, then you can run the
    ```make```
    command that will automatically produce the .exe for you. **But be careful**: it is very likely that you will need to compile your program with very specific compilator and/or linker flags (like those needed for openmp, for instance), so make sure to add all the flags you need. After the compilation the executable "fortran_tester" is produced.
