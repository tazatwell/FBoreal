REM !\bin\bash
REM Run something (anything) using the simulator just to verity that it runs
REM ASSUMES the FBP ".so" file is in this directory
REM IMPORTANT: the Dogrib data will use too much memory on most computers

SET DDIR=%C:\Users\dlm\FBoreal\Instances\9CellsHom\%

SET TDIR=%C:\Users\dlm\FBoreal\FBoreal\%

python %TDIR%SimApril2018.py --input-instance-folder=%DDIR% --ignitions --nsims=1 --output-grid --plot --combine --ROS-Threshold=0.1 --output-folder=outdir --weather=rows





