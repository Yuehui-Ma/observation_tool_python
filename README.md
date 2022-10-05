# observation_tool_python

## A script can be used to make observational plans for the Purple Mountain Observatory (PMO) 13.7m mm-wavelength telescope.

### Input: A txt file recording the central coordinates of the candidate off positions (in decimal l, b format)

The PMO-13.7m telescope uses the position-switch On-The-Fly mode for mapping observations. The sky area is usually divided into units of 0.5 arcdeg $\times$ 0.5 arcdeg. The off positions are those areas, 10 $\times$ 10 arcmin^2 in size, in the vicinity of the source units without emissions above a given sensitivity. Users can use the obs_tool.py script to prepare for the observations. The tool is provided by an input 'txt' file, which records the candidate off positions, and can calculate the observable areas using these off positions within a specific LST time interval. The Default observable condition is that the source altitude should be higher than 30 arcdeg and lower than 80 arcdeg, the altitude difference between the source and off position should be less than 1 arcdeg, and their azimuthal difference should be less than 5 arcdeg. 

### Outputs: 
- A '.pdf' figure marks the locations of the off areas and the observable sky areas, with each color of marker corresponding to an off position. (G120_obs.pdf is an example)
- Several '.cat' files (depending on the numbers of input off positions) record the source names of the observable areas and the observable LST time intervals.  (F122.750+08.225.cat is an example)
