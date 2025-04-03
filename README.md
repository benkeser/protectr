# Code repository for PROTECT

## Getting started on RSPH HPC

### Initial setup

Follow these steps to set up the repository and environment on RSPH HPC.

1. **Clone the repository** into your home directory:

   ```bash
   git clone git@github.com:benkeser/protectr.git
   ```

2. **Create a personal R library and install required packages**:

   ```bash
   # Create a directory for R libraries
   mkdir ~/Rlibs

   # Start an interactive session
   srun --pty --partition=interactive-cpu --nodes=1 --ntasks-per-node=1 --mem-per-cpu=8G --time=02:00:00 bash

   # Load R module
   module load R/4.4.0

   # Open R console
   R

   # Set library path and install packages
   .libPaths("~/Rlibs")

   # Choose a USA mirror when prompted
   install.packages(c("here", "tidyverse", "fastverse", "future.apply", "progressr", "config", "assertthat"))
   ```

3. **Transfer data files** to the `protectr/data` folder using `scp` or an SFTP client (e.g., CyberDuck).

4. **Create directories for results and scratch files in the project space**:

```bash
cd /projects/dbenkes
mkdir seth

cd seth
mkdir protectr

cd protectr

mkdir boot_res
mkdir scratch
```

5. **Update output paths**:

- Change output paths in `run_simulation.sh` and `run_bootstrap.sh` to the scratch folder created in step 4. This can be done using vim or nano text editor on HPC directly, or make changes locally then push to GitHub. To use vim:
 
 ```bash
cd ~/protectr
vi run_simulation.sh

# Use keypad to move cursor down to output line
# Click 'i' to enter insert mode
# Swap 'allison' for 'seth'
# Click 'ESC'
# Click ':wq' for write (save) and quit

# Repeat for run_bootstrap.sh

```


6. **Update result paths**:

- Change path to save boostrap results in `run_bootstrap.R` to the boot_res folder created in step 4. This can be done using vim or nano text editor on HPC directly, or make changes locally then push to GitHub. To use vim:

```bash
cd ~/protectr
vi run_bootstrap.R

# Use keypad to move cursor down to final line that saves the RDS file
# Click 'i' to enter insert mode
# Swap 'allison' for 'seth'
# Click 'ESC'
# Click ':wq' for write (save) and quit

```

- Change path to load the boostrap results and create a bootstrap confidence interval in `run_bootstrap_ci.R` to the boot_res folder created in step 4. This can be done using vim or nano text editor on HPC directly, or make changes locally then push to GitHub. To use vim:

```bash
cd ~/protectr
vi run_bootstrap.R

# Use keypad to move cursor down to final line that saves the RDS file
# Click 'i' to enter insert mode
# Swap 'allison' for 'seth'
# Click 'ESC'
# Click ':wq' for write (save) and quit

```

### Creating weekly records data and obtaining initial point estimates

Submit a job using `run_simulation.sh` to generate weekly records and obtain point estimates. For example, to run a job for Haiti on David's partition, you would submit a job as follows:

```bash
cd ~/protectr
./run_simulation.sh benkeser haiti
```


- Check job status with:
  ```bash
  squeue -u <your_user_id>
  ```
- It takes \~2 hours to create the weekly records data (varies by site). If weekly records data for a givens site already exists, it will load it directly/skip this step.
- Fitting propensity models and MSMs runs more quickly (\~20 minutes)
- Outputs/errors are logged in `/projects/dbenkes/seth/protectr/scratch/<site>_<jobid>.out`.
- Results are saved in `~/protectr/results/<country_name>`.

### Running boostrap and obtaining confidence intervals

Submit a bootstrap job with `run_bootstrap.sh` to run `run_bootstrap.R` n_boot number of times and save results in the project space. For example, to run a job for Haiti on David's partition for 1000 seeds/bootstrap replicates, you would submit a job as follows:

```bash
./run_bootstrap.sh benkeser haiti 1000
```

Once all bootstrap replicates are complete, aggregate results:

```bash
srun --pty --partition=interactive-cpu --nodes=1 --ntasks-per-node=1 --mem-per-cpu=8G --time=02:00:00 bash

module load R/4.4.0

Rscript run_bootstrap_ci.R haiti
```

