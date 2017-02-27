**These files are the python source code for the Dengue Model. There are 12 total *.py files:**

- dengue_run.py 
	- This is the main wrapper script that runs the main model simulation (dengue_model.py), parameter estimation (param_est.py), and profile likelihood (profile.py and residual_sigma.py). 
	- dengue_model.py		
	- param_est.py	
	- profile.py
	- residual_sigma.py
- plot_proflike_general.py
	- Visualizes the profile likelihood outcome from dengue_run.py
- fim_test.py
	- This script runs the Fisher information matrix (fim.py) for the main model simulation
	- fim.py
- pair_surface.py 
	- This script calculates the likelihood surface when iteratively changing pairs of parameters
- heatmap_loop_pair.py
	- Visualizes the likelihood surface from pair_surface.py
- sim.py
	- Runs a dengue intervention (e.g., decreasing the number of mosquito larvae) using dengue_model_intervention.py
- dengue_model_intervention.py	




> Written with [StackEdit](https://stackedit.io/).