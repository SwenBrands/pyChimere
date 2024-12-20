# pyChimere

A Python package to run, develop and validate the chemical transport
model CHIMERE, version 2017, on a High Performance Cluster (CESGA's
Finisterrae 3 in this case). This includes the proper processing of 
chemical species concentrations at the lateral boundary conditions
provided by C-IFS, visualiztion software, and scripts to merge
anthropogenic emssions from the EMEP inventory with point-wise PRTR
emission sources gathered by the Galician Government (Xunta de Galicia).

Author: Swen Brands, brandssf@ifca.unican.es or swen.brands@gmail.com

Reference article:  https://doi.org/10.5194/gmd-13-3947-2020

Chimere home page: https://www.lmd.polytechnique.fr/chimere/


# The bash scripts in this directory have the following funtions #######

1. lanza_chimere.sh: sends the op_nueva.sh script from the frontal node
to the working nodes

2. variables_default.sh: the configuration file called by lanza_chimere.sh

3. op_nueva.sh: the main script containing the model pre-processing steps,
the model execution itself and model postprocessing tasks

4. runme.sh: run a large series of experiments for different init dates



# The directories cover the following tasks ############################

1. copernicus: scripts to post-process oprational C-IFS predictions so
that they can be fed to CHIMERE at its lateral boundaries

2. emissions: scripts to merge EMEP inventory with local emissions from
the PRTR archive hosted by the Galician Government (MeteoGalicia)

3. figures: scripts to map model results, model parameters and emission
output of the emisurf prgramm provided with chimere. 

4. validation: scripts to validate the model output against in-situ
observations from the Galician "Red de Calidad del Aire"


Credits
-------
This ongoing research line is being funded by Xunta de Galicia.

<img align="right" width="200" src="https://www.xunta.gal/ficheiros/identidade-corporativa/2021/simbolos/simbolo-positivo.svg">
