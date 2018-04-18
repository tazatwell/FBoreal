REM !\bin\bash
REM Run something (anything) using the simulator just to verity that it runs
REM ASSUMES the FBP ".so" file is in this directory
REM IMPORTANT: the Dogrib data will use too much memory on most computers

SET DDIR=%C:\Users\dlm\FBoreal\Instances\9CellsHom\%
REM Note that I inserted the \ after 9CellHOM because if it was not there input data file, it kept looking for 
REM a file named \9CellsHomWeather.csv instead of \9CellsHom\Weather.csv becuase of the way the %DDIR% environment variable
REM was referenced in the python command line
REM This is consistet with the way the source file (SimOct2017) is referenced in the python command line

SET TDIR=%C:\Users\dlm\FBoreal\FBoreal\%

python %TDIR%SimApril2018.py --input-instance-folder=%DDIR% --ignitions --nsims=1 --output-grid --weather-file --plot --combine --ROS-Start-Threshold=1e-3 --Mean-Burnout-Periods=2 --output-folder=outdir





