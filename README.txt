## Installation of PlanarFold (on Linux)

1. Install Miniconda（if not installed）：
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  $ bash Miniconda3-latest-Linux-x86_64.sh

2. Create and activate conda environment：(this step might take a while)
  $ conda env create -f environment.yml
  $ conda activate planarfold

3. Compile the software：
  $ cd PlanarFold/src
  $ chmod +x build.sh
  $ ./build.sh

4. Execute the software：
  $ cd PlanarFold/tools
  $ ./PlanarFold -h

5. Set up environment:
  For bash, add the following lines into the "~/.bashrc" file
  ''''
  # PlanarFold
  PATH=$PATH:/xxxx/PlanarFold/tools
  export PYTHONPATH=$PYTHONPATH:/xxxx/PlanarFold/src
  ''''
  The 'xxxx' is the directory of PlanarFold

  For csh, add the following lines into the "~/.cshcr" file
  ''''
  # PlanarFold
  setenv PYTHONPATH /xxxx/PlanarFold/src
  setenv PATH ${PATH}:/xxxx/drfold_final_TR/PlanarFold/tools
  ''''
  The 'xxxx' is the directory of PlanarFold
