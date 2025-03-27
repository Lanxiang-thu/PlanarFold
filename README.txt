## Installation of PlanarFold (on Linux)

1. Install Miniconda（if not installed）：
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  $ bash Miniconda3-latest-Linux-x86_64.sh

2. Create and activate conda environment：(this step might take a while)
  $ conda env create -f environment.yml
  $ conda activate my_python3.9

3. Compile the software：
  $ cd PlanarFold/src
  $ chmod +x build.sh
  $ ./build.sh

4. Execute the software：
  $ cd PlanarFold/tools
  $ ./PlanarFold -h
