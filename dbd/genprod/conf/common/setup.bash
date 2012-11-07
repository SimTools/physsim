##
##  Setup environments to run root/jsf/Geant4
##  Environment parameter for LC studies
##
#######################################################
##
##  Setup for JLCLOGIN2
##
#######################################################

export LC_RELEASE=1.46
# export LC_RELEASE_DIR=/proj/soft/Release/$LC_RELEASE
# export LC_RELEASE_DIR=/home/ilc/miyamoto/soft/MyRelease/110421
# export LC_RELEASE_DIR=/home/ilc/miyamoto/soft/MyRelease/110829
export LC_RELEASE_DIR=/home/ilc/miyamoto/soft/MyRelease/121030
#
utildir=/group/ilc/soft/utils64
. ${utildir}/gf44.setup
export ROOTSYS=${utildir}/root/5.28.00h
export G4INSTALL=${utildir}/g4/geant4.9.5.p01
# export LCIO=${utildir}/../ilcsoft/x86_64-sl5/v01-09/lcio/v01-51-02
export LCIO=${utildir}/lcio/v02-01-02
export JDK_HOME=/usr/java/jdk1.7.0_01
export CERN_ROOT=${utildir}/cernlib/2005

export CLHEP_BASE_DIR=${utildir}/CLHEP/2.1.1.0
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export STDHEPDIR=${utildir}/stdhep/5.06.01

export LCBASEDIR=$LC_RELEASE_DIR/lcbase
export LEDAROOT=$LC_RELEASE_DIR/Leda
export LCLIBROOT=$LC_RELEASE_DIR/lclib
export JSFROOT=$LC_RELEASE_DIR/jsf
export KFLIBROOT=$LC_RELEASE_DIR/physsim

export JUPITERROOT=$LC_RELEASE_DIR/Jupiter
export SATELLITESROOT=$LC_RELEASE_DIR/Satellites
export URANUSROOT=$LC_RELEASE_DIR/Uranus
export SOSYMLINK=true


# Define following environment parameters to build JSF with Whizard 
# for DBD studies
export WHIZDIR=/group/ilc/soft/utils64/whizard/whizard-1.95
# export utiltau=/home/ilc/miyamoto/soft/physsim2
export utils=/nfs/g/ilc/soft/utils64
export TAUOLADIR=${utildir}/tauola/tauola_desy/TAUOLA/tauola
export PHOTOSDIR=${utildir}/tauola/tauola_desy/TAUOLA/photos
export PYTHIADIR=${utildir}/pythia/6.4.22

##### Geatn4 setup ##############
# . $G4INSTALL/env.sh > /dev/null
######### end of Geant4 setup ##################################

## Set command path
LCPATH=$LCBASEDIR/bin:$JSFROOT/bin:$ROOTSYS/bin:$CERN_ROOT/bin:$LCLIBROOT/bin:$CLHEP_BASE_DIR/bin
# J4PATH=$JUPITERROOT/bin/$G4SYSTEM:$G4INSTALL/bin/$G4SYSTEM
LCIOPATH=$LCIO/tools:$LCIO/bin:$JDK_HOME/bin
export PATH=.:~/bin:$LCPATH:$LCIOPATH:$PATH
unset LCPATH
unset J4PATH
unset LCIOPATH

## Set LD Library Path
export LD_LIBRARY_PATH=$JSFROOT/lib:$ROOTSYS/lib:$LEDAROOT/lib:$CLHEP_BASE_DIR/lib:${KFLIBROOT}/lib:${LD_LIBRARY_PATH}
export IMAKEINCLUDE="-I$LCBASEDIR -I$KFLIBROOT -I$LCLIBROOT"

