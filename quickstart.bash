#b!\bin\bash
# Run something (anything) using the simulator just to verity that it runs
# ASSUMES the FBP ".so" file is in this directory
# IMPORTANT: the Dogrib data will use too much memory on most computers

BASEDIR="/home/woodruff/Documents/Research/FBoreal"
DDIR="$BASEDIR/Instances/9CellsHom/"
TORUN="$BASEDIR/SimApril2018.py"

python $TORUN --input-instance-folder=$DDIR --ignitions --nsims=1 --output-grid --plot --combine --ROS-Threshold=0.1 --output-folder=outdir --weather=rows

##python SimJan2018.py --input-instance-folder=$DDIR --output-folder=simoutput --ignitions --nsims=1 --output-grid --weather-file --plot --combine --ROS-Start-Threshold=1e-3 --spotting-parameter-data-file-name=SpotSample.json
