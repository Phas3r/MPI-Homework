@echo off
SET path_to_mpi="C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec.exe"
SET path_to_app="C:\Users\Sergey\Documents\Visual Studio 2015\Projects\ConsoleApplication1\x64\Release\ConsoleApplication1.exe"

SET upper_limit=1.0
for /l %%m in (1, 1, 4) do (
    for /l %%x in (1, 1, 20) do (
        for /l %%i in (1,1,50) do (
            %path_to_mpi% -n %%x %path_to_app% %%m %upper_limit% >> results%%m.txt
        )
    )
)