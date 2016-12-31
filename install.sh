#!/usr/bin/env bash
# Script to setup camoco in an anaconda environment
# Written by Joe Jeffers and Rob Schaefer
# Email: jeffe174@umn.edu, schae234@umn.edu


function usage(){
cat <<EOF
    Usage: $0 [flags]

    Flags
    -----
    -h | --help
        print help message
    -b | --base
        Base installation directory for camoco (default: ~/.camoco).
EOF
exit 0
}

# Configurable variables
BASE=$HOME/.camoco

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

function red(){
    printf "${RED}$1${NC}\n"
}
function green(){
    printf "${GREEN}$1${NC}\n"
}

while [[ $# > 0 ]]
do
key="$1"
case $key in 
    -h|--help)
    usage
    shift
    ;;
    -b|--base)
    BASE=$2
    shift
    ;;
    *)  
        #unknown options
    ;;
esac
shift
done

NAME="camoco"

#===================================================
#----------Check for internet connection------------
#===================================================
wget -q --tries=10 --timeout=20 --spider http://github.com
if [[ $? -ne 0 ]]; then
        echo "Check your internet connection and try again."
		exit
fi

export CWD=$(pwd)
#===================================================
#----------Setup the build Environment--------------
#===================================================
echo "Setting up the build environment"
source $HOME/.bashrc
mkdir -p $BASE
mkdir -p $BASE/conda
mkdir -p $BASE/bin 
mkdir -p $BASE/lib 

cd $BASE

export LD_LIBRARY_PATH=$BASE/lib:$LD_LIBRARY_PATH
export PATH=$PATH:$BASE/bin:$BASE/conda/bin

#===================================================
#----------------Install conda ---------------------
#===================================================
# if ! hash conda 2>/dev/null
if [ ! -e $BASE/conda/bin/conda ]
then
    cd $BASE
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -f -p $BASE/conda
	rm -f miniconda.sh
else
    green "Conda Already Installed" 
fi

#===================================================
#----------------Install git lfs--------------------
#===================================================
if ! hash git-lfs 2>/dev/null ; then
    echo "Installing git-lfs for the large files in the repo"
    cd $BASE
    wget https://github.com/github/git-lfs/releases/download/v1.1.0/git-lfs-linux-amd64-1.1.0.tar.gz
    tar xzf git-lfs-linux-amd64-1.1.0.tar.gz
    rm -rf git-lfs-linux-amd64-1.1.0.tar.gz
    cd git-lfs-1.1.0/
    rm -rf $BASE/bin/git-lfs*
    mv git-lfs $BASE/bin/
    cd $BASE
    rm -rf git-lfs*
    git lfs install
    git lfs uninstall
else
    green "Git-lfs already installed"
fi

#===================================================
#-----Install libraries to Ensure Smoothness--------
#===================================================
if [ ! -e $BASE/bin/icu-config ]
then
    echo "Installing a local version of the icu library"
    cd $BASE
    wget http://archive.ubuntu.com/ubuntu/pool/main/i/icu/icu_4.8.1.1.orig.tar.gz
    echo "Extracting the tarball"
    tar xzf icu_4.8.1.1.orig.tar.gz
    rm -rf icu_4.8.1.1.orig.tar.gz
    cd icu/source
    ./configure --prefix=$BASE
    make
    make install
    cd $BASE
    rm -rf icu/
else
    green "icu already installed"
fi
if [ ! -e $BASE/lib/libqhull.so ]
then
    echo "Installing a local version of the qhull library"
    cd $BASE
    wget http://www.qhull.org/download/qhull-2012.1-src.tgz
    echo "Extracting the tarball"
    tar xzf qhull-2012.1-src.tgz
    rm -rf qhull-2012.1-src.tgz
    cd qhull-2012.1
    make
    mv bin/* $BASE/bin/
    mv lib/* $BASE/lib/
    ln -s $BASE/lib/libqhull.so $BASE/lib/libqhull.so.5
    cd $BASE
    rm -rf qhull-2012.1/
else
    green 'qhull already installed'
fi

#===================================================
#------------------Install mcl----------------------
#===================================================
if [ ! -e $BASE/bin/mcl ]
then
    echo "Installing mcl locally"
    cd $BASE
    wget http://micans.org/mcl/src/mcl-latest.tar.gz
    echo "Extracting the taball"
    tar xzvf mcl-latest.tar.gz
    rm -rf mcl-latest.tar.gz
    cd $(find . -name 'mcl-*' -type d | head -n 1)
    ./configure --prefix=$BASE
    make
    make install
    cd $BASE
    rm -rf $(find . -name 'mcl-*' -type d | head -n 1)
else
    green 'mcl already installed'
fi

#===================================================
#----------Build the Conda Environment--------------
#===================================================
if [ ! -d $BASE/conda/envs/camoco ]
then
    echo "Making the conda virtual environment named $NAME in $BASE"
    conda remove -y --name $NAME --all
    conda config --add envs_dirs $BASE/conda/envs
    conda config --append channels conda-forge
    conda create -y -n $NAME python=3 setuptools pip cython numpy scipy pandas \
        matplotlib feather-format nose six pyyaml yaml pyparsing python-dateutil \
        pytz numexpr patsy statsmodels networkx mpmath termcolor scikit-learn \
        ipython ipdb pytest-cov flask gunicorn
else
    green 'conda already installed'
fi

#===================================================
#----------Activate the Conda Environment-----------
#===================================================
green "activating $NAME"
source $BASE/conda/bin/activate $NAME
green 'checking python'
which python

#==================================================
#----------Take care of some pip packages ---------
#==================================================
python -c 'import powerlaw'
if [ $? -eq 1  ]
then
    pip install powerlaw
else
    green "powerlaw installed"
fi

#===================================================
#-----------------Install apsw----------------------
#===================================================
python -c 'import apsw'
if [ $? -eq 1 ]
then
    echo "Installing apsw into the conda environment"
    cd $BASE
    rm -rf apsw
    git clone https://github.com/rogerbinns/apsw.git
    cd apsw
    python setup.py fetch --missing-checksum-ok --all build --enable-all-extensions install
    cd $BASE
    rm -rf apsw
else
    green "apsw installed"
fi


#==================================================
#---------------Update Setuptools------------------
#==================================================
pip install setuptools --upgrade

#==================================================
#-----------------Install Camoco-------------------
#=================================================
green "Installing Camoco"
cd $CWD
python setup.py install
python -c 'import camoco'
if [ $? -eq 1 ]
then
    red 'Camoco failed to install!'
    exit 1
else
    green 'Camoco installed!'
fi
source deactivate 

#===================================================
#------------Update the bashrc----------------------
#===================================================
if [ $(grep $BASE/conda/bin ~/.bashrc | wc -l) -eq 0 ]
then 
    red '-----------------------------------------------'
    red "Update your $HOME/.bashrc:"
    echo "export LD_LIBRARY_PATH=$BASE/lib:\$LD_LIBRARY_PATH"
    echo "export PATH=$BASE/bin:$BASE/conda/bin:\$PATH"
    red '-----------------------------------------------'
fi

#===================================================
#-------------Use Instructions----------------------
#===================================================
green '==============================================='
echo "All done, to use camoco, first restart the shell."
echo "Then type: source activate $NAME"
echo 'Start with the Camoco manual page: camoco --help'
echo "When finished, type: source deactivate"
echo ""
echo ""
echo "Please report all bugs to: https://github.com/schae234/Camoco/issues"
green '==============================================='
