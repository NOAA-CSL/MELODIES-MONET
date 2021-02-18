These directories are model specific
	cmaq
	rapchem
	wrfchem

We are combining these to create one code for cmaq, fv3, wrfchem, and rapchem. 
Also make it easy to add other models and mechanisms too
	unified 

List of components needed in a namelist or to add to bash scripts
	-radius of influence (dependent on model resolution)
	-model
	-output name
	-mapping table

Progress
	Becky adapted some of Patrick's scripts to be more generalizable with the above components. 
	Would other approaches like a namelist be better?
	Becky wants adaptability for reading in multiple files to average over month time scale, which is available with current bash scripts
