#!/bin/bash
#SBATCH --job-name=compile      # nom du job
#SBATCH --ntasks=1
#SBATCH --time=00:20:00            # Temps d’exécution maximum demande (HH:MM:SS)
#SBATCH --output=out_compile%j.out  # Nom du fichier de sortie
#SBATCH --error=err_compile%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --partition=compil
#SBATCH -A rkm@cpu


# print all commands
set -x

# cleanup
module purge

#==============================
# STEP 1: COMPILE GLFW
# -> make sure DBUILD_SHARED_LIBS=ON
# -> files https://github.com/glfw/glfw.git

if [ 0 ]; then

	module load cmake/3.31.4
	module load gcc/14.2.0

	cd $HOME
	mkdir lib; cd lib;
	git clone https://github.com/glfw/glfw.git
	cd glfw
	rm -r build # to reset any previous cmake variable declaration
	kdir build

	C=gcc cmake -S . -B build \
	 -DCMAKE_INSTALL_PREFIX=$HOME/lib/glfw/install \
	 -DBUILD_SHARED_LIBS=ON \
	 -DGLFW_BUILD_EXAMPLES=OFF \
	 -DGLFW_BUILD_TESTS=OFF \
	 -DGLFW_BUILD_DOCS=OFF \
	 -DGLFW_BUILD_WAYLAND=OFF \
	 -DGLFW_BUILD_X11=OFF

	make --build build
	make --install build

	module purge
	echo 'Done compiling GLFW'
fi
#==============================
# STEP 2a: COMPILE QCC

module load gcc/14.2.0
module load cuda/12.6.3

export BASILISK=$HOME/basilisk/src
cd $BASILISK

# qcc, bview
ln -sf config.gcc config # fixme: use proper JeanZay config
make qcc

echo 'Done compiling qcc'
#=============================
# STEP 2b: COMPILE GPU BACKENDS

if [ 0 ]; then
	# gpu (opengl)
	# I need glfw.h 
	cd $BASILISK/grid/gpu
	export LIBGLFW=$HOME/lib/glfw/install # STEP 1
	INCLUDE="${LIBGLFW}/include"
	CFLAGS="-I${INCLUDE}"
	CFLAGS=$CFLAGS make
	echo 'Done compiling opengl backend'
fi

# gpu (cuda)
cd $BASILISK/grid/cuda
CFLAGS=$CFLAGS make libbuda.a
echo 'Done compiling cuda backend'


