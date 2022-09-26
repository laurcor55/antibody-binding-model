# Brownian Diffusion Binding Model
This repository contains all the code associated with the Brownian Diffusion binding model. This file is currently specific to MacOS.

## Installation
1. **Install Python3** - To check if python3 is installed, open Terminal app. Type `python3 --version` in terminal. If installed, it will say the python version. For example, `Python 3.9.2`. If not installed, install from [Python website](https://www.python.org/).

2. **Install Pip** - Make sure Pip is installed by running `python3 -m ensurepip --upgrade` in the Terminal app.

3. **Download Model Folder** - In repository, click green "Code" dropdown menu at top right. Select "Download Zip" to download folder. Move and unzip model folder where desired.

4. **Open Model Folder in Terminal** - Drag and drop model folder into Terminal application. This should navigate to the model folder.

5. **Install Requirements** - In the Terminal app, run `pip install -r requirements.txt`.

## Usage
1. **Open Model Folder in Terminal** - Drag and drop model folder into Terminal application. This should navigate to the model folder.

2. **Run Application** - Run command `python3 run.py`. 

3. **Configure Parameters** - Change parameters with "Ligand Input", "Substrate Input" and "Reaction Count" tabs on the left panel, towards the center. To see setup, press "Preview Reaction" and click left plot. Red sphere is the ligand. Black spheres are the substrates.

4. **Start Reaction** - Press "Start Reaction" button to run reaction simulations. If needed, you can stop the simulation by pressing "Stop Reaction".

5. **View Output** - Model output seen on the right panel. 

6. **Close App** - Press top X to exit app. 