#==============================================================================
#     Download and compile recent Petsc, then source its env. vars to .bashrc
#     You must have OpenFOAM vars sourced and sudo privillege
#==============================================================================

# Install dependencies
sudo apt update
sudo apt-get install libatlas-base-dev libblas-dev \
    liblapack-dev flex bison git make cmake gfortran

set -vx

# Choose latest PetSc version
# Ubuntu repos have 3.7.7
PetscV="3.14.0"

# Install as OpenFOAM's ThirdParty software
PetscIntallDir=$WM_PROJECT_USER_DIR/ThirdParty
mkdir -p $PetscIntallDir
cd $PetscIntallDir

# MPI Paths
mpiDir=$MPI_ARCH_PATH
mpiFort=$(which mpifort) 
mpiCC=$(which mpicc)
mpiCXX=$(which mpicxx)

## Download Petsc (Comment out this part if installing from local disk)
[ -z "$WM_PROJECT_DIR" ] && \
    { echo "\nERROR: Please source OpenFOAM's /etc/bashrc\n"; exit 1; }
echo "\nDownloading Petsc v$PetscV \n"
wget -c "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PetscV.tar.gz"

# Untar the file
tar -xvzf petsc-$PetscV.tar.gz

## Configure petsc
cd petsc-$PetscV 
./configure --with-mpi-dir=$mpiDir --with-debugging=0 \
    --with-precision=double --with-shared-libraries=1 --with-scalar-type=real \
    --with-fc=$mpiFort --with-cc=$mpiCC --with-cxx=$mpiCXX \
    CXX_LINKER_FLAGS=-Wl,--no-as-needed \
    CFLAGS="-g -O2 -fPIC -fstack-protector-strong -Wformat -Werror=format-security" \
    CXXFLAGS="-g -O2 -fPIC -fstack-protector-strong -Wformat -Werror=format-security" \
    FCFLAGS="-g -O2 -fstack-protector-strong" \
    FFLAGS="-g -O2 -fstack-protector-strong" \
    CPPFLAGS="-Wdate-time -D_FORTIFY_SOURCE=2" \
    LDFLAGS="-Wl,-Bsymbolic-functions -Wl,-z,relro" MAKEFLAGS=w
#    --download-hypre --download-parmetis --download-metis \
#    --download-ptscotch --download-mumps --download-scalapack \

echo "\nCompiling Petsc ..."

parch=$(tail -200 configure.log | grep 'PETSC_ARCH:' | awk '{print $2}')
pdir=$(tail -200 configure.log | grep 'PETSC_DIR:' | awk '{print $2}')

## Compile petsc
make PETSC_DIR=$pdir PETSC_ARCH=$parch all

## Source Petsc env. vars.
echo "export PETSC_ARCH=$parch" >> $HOME/.bashrc
echo "export PETSC_DIR=$pdir" >> $HOME/.bashrc
echo "export LD_LIBRARY_PATH=\$PETSC_DIR/\$PETSC_ARCH/lib:\$LD_LIBRARY_PATH" >> $HOME/.bashrc
echo "\nOK."
