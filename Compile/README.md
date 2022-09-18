cd to `/Compile` folder in your HPC. Upload `Borg-1.9` folder to `Compile` folder. You can download `Borg-1.9` from http://borgmoea.org/#contact and you need to save the Borg folder under this directory

Also upload the `UHMOEA_ms.cpp` and makefile files to this folder. Then run make. This will create `borgms.o`, `mt19937ar.o`, `UrbanHeat.o`, and `UrbanHeatMPI` that are used for optimization.

Once all files are prepared, run `make`.

You might need to type in `chmod 755 UrbanHeatMPI` to get permission for using MPI on your machine.

Now, proceed to the `/Optimization` folder and consult the `README.md` in that directorty.
