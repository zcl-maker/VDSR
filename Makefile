                        #
##########################################################################################################

#========================================== Program name ================================================#
PN = vdsr_sjtu.e
#========================================================================================================#

#===================================== Postprocessor program name =======================================#
PPN = vdsr.e
#========================================================================================================#


#=========================================== Define V_MPI ===============================================#
DEF_V_MPI=-DV_MPI
#========================================================================================================#

#=========================================== Define Parallel HDF5 =======================================#
#========================================================================================================#

#======================================== Define SYNCHROTRON ============================================#
DEF_SYNCHROTRON=-DSYNCHROTRON
#========================================================================================================#


#========================================= mpich directory ==============================================#
#MPI_DIR=/opt/mpich/ch_p4
#MPI_DIR=/opt/mpt/3.5.0/xt/mpich2-pgi
MPI_DIR=/cluster/soft/apps/mvapich/1.0.1/gcc.pgf90
#========================================================================================================#

#============================== Message Passing Interface include files =================================#
#MPI_INCLUDE = -I/.netmount/storage_giga/$(MPICH)/ch_p4/include/
#MPI_INCLUDE = -I/.netmount/storage_giga/$(MPICH)/ch_p4mpd/include/
MPI_INCLUDE = -I$(MPI_DIR)/include
#========================================================================================================#

#========================================== MPI lib files ===============================================#
#MPI_LIB_DIR = -L/.netmount/storage_giga/$(MPICH)/ch_p4/lib/
#MPI_LIB_DIR = -L/.netmount/storage_giga/$(MPICH)/ch_p4mpd/lib/
MPI_LIB_DIR = -L$(MPI_DIR)/lib/

MPI_LIB = -lmpich

#========================================= HDF5 directory ==============================================#
HDF5=/cluster/soft/hdf5/1.6.10/gcc.ifort
#HDF5=/opt/hdf5
#HDF5=/opt/cray/hdf5/1.8.3.0/hdf5-pgi
#HDF5=/usr/common/usg/hdf5/1.6.5/parallel
# the next one is for Franklin
#HDF5=/opt/cray/hdf5-parallel/1.8.3.1/hdf5-parallel-pgi
#HDF5=$(HDF5_DIR)
# the next one is for Hopper
#HDF5=/opt/cray/hdf5-parallel/1.8.5.0/hdf5-parallel-pgi
#========================================================================================================#

#============================== Message Passing Interface include files =================================#
HDF5_INCLUDE = -I$(HDF5)/include/
#========================================================================================================#

#========================================== MPI lib files ===============================================#
HDF5_LIB_DIR = -L$(HDF5)/lib/
HDF5_LIB = -lhdf5 -lz
#========================================================================================================#

#======================================== PulseProfiles lib =============================================#
PULSE_PROFILE_LIB_DIR = PulseProfile
PULSE_PROFILE_LIB = -lPulseProfiles
#========================================================================================================#

#======================================== PulseProfiles lib =============================================#
SCALES_LIB_DIR = Scales
SCALES_LIB = -lScales
#========================================================================================================#

#=============================== GNU Scientific library ============= ===================================#
#GSL_DIR = /opt/gsl/
GSL_DIR= /usr/common/usg/gsl/gsl-1.10
#=============================== GNU Scientific library include files ===================================#
#GSL_INCLUDE = -I/.netmount/storage_vlpl/usr/include/gsl/
#GSL_INCLUDE = -I/.netmount/storage_giga/usr/include/gsl/ -I/.netmount/storage_giga/usr/include/
GSL_INCLUDE = -I$(GSL_DIR)/include
#========================================================================================================#

#================================= GNU Scientific library lib files =====================================#
#GSL_LIB_DIR = -L/.netmount/storage_vlpl/usr/lib/
#GSL_LIB_DIR = -L/.netmount/storage_giga/usr/lib/
GSL_LIB_DIR = -L$(GSL_DIR)/lib/
GSL_LIB = -lgsl -lgslcblas
#========================================================================================================#

#======================================= KTrace include files ===========================================#
#KTRACE_INCLUDE = -I/usr/include/kde/
#========================================================================================================#

#======================================== OLDVERSION Define =============================================#
# It is define for compatibility with old version of saving files.
#DEFINE_OLDVERSION=-DOLDVERSION
#========================================================================================================#

##########################################################################################################
#============================================== Make ====================================================#
M_MAKE = $(MAKE) -j4
#M_MAKE = $(MAKE)
#========================================================================================================#

#============================================ Compiler ==================================================#
#CC = xlC
#mpCCFlags = -qrtti
#CC = g++
#VLPLCC = icpc # for debug version with Intel compiler
VLPLCC=mpiCC
#the up CC=/opt/cray/xt-asyncpe/3.3/bin/CC
#CC = icc # Intel compiler
#FS = -ffloat-store
#NO-DEPRECATED = -Wno-deprecated
#WARNING = -W -Wall
#========================================================================================================#

#====================================== Intel Compiler flags ============================================#
# If you use Intel compiler please specify additional flags
#ICC_FLAGS=-static-libgcc -mtune=pentium4 -xW -ip -parallel $(WARNING)  //canceled by MChen 091110
GCC_FLAGS= 
VLPLCC_FLAGS=-DMPICH_IGNORE_CXX_SEEK
#========================================================================================================#

#========================================= Compiler flags ===============================================#
Optimized = -O2
Debug = -O0
CFLAGS = $(mpCCFlags) $(OO) $(GG) $(VLPLCC_FLAGS) -c $(MPI_INCLUDE) $(GSL_INCLUDE) $(HDF5_INCLUDE) \
	 $(KTRACE_INCLUDE) $(FS) $(DEFINES) \
         $(MPI) $(NO-DEPRECATED) $(WARNING)
STATIC = -static
#========================================================================================================#

#=========================================== Link flags =================================================#
LINK = $(VLPLCC) $(VLPLCC_FLAGS) -o # if you use icc then specify "ICC_FLAGS"
LINK_FLAGS = $(LF)
#========================================================================================================#

#======================================= vlpl3d project libs ============================================#
PROJECT_LIB_DIR = lib
PROJECT_LIB_DIR_FLAG = -Llib ${MPI_LIB_DIR} ${HDF5_LIB_DIR} ${GSL_LIB_DIR}
PROJECT_LIB_DIR_FLAGS =${PULSE_PROFILE_LIB_DIR} ${SCALES_LIB_DIR} ${MPI_LIB_DIR} ${HDF5_LIB_DIR} ${GSL_LIB_DIR}
#PROJECT_LIBS_FLAGS =${PULSE_PROFILE_LIB} ${SCALES_LIB} ${GSL_LIB} ${HDF5_LIB}
PROJECT_LIBS_FLAGS =${HDF5_LIB}
#========================================================================================================#

#========================================== lib utilities ===============================================#
AR_LIB = /usr/bin/ar
AR_FLAGS = rv
RAN_LIB = /usr/bin/ranlib
#========================================================================================================#

#===================================== vlpl3d project libs define =======================================#
P_LIBS = vlpl3d_lib
#========================================================================================================#

#======================================== vlpl3d project files ==========================================#

#======================================== vlpl3d project files ==========================================#
objects =  buffers.o Cgs.o Detector.o domain.o myhdfshell.o SpecialFunc.o\
		namelist.o Particle.o Pixel.o vdsr.o Beam.o ExternalFields.o Controls.o Photon.o
#========================================================================================================#

#====================================== vlpl3d project objects ==========================================#
objects_vlpl3d = $(objects)
objects_vlpl3d_link = $(objects)
#========================================================================================================#
##########################################################################################################

#========================================== vlpl3d project ==============================================#
vl : vlpl3d strip

vlpl3d :
	$(M_MAKE) OO=$(Optimized) LF="${PROJECT_LIBS_FLAGS} ${MPI_LIB} ${PROJECT_LIB_DIR_FLAG}" $(PN) \
   MPI=$(DEF_V_MPI) DEFINES="${DEFINE4ALL}"

$(PN) : $(objects_vlpl3d)
	${LINK} $(PN) $(objects_vlpl3d_link) $(LINK_FLAGS)
#========================================================================================================#

#====================================== vlpl3d project DEBUG ============================================#
vld : makevld

makevld :
	$(M_MAKE) OO=$(Debug) GG=-g LF="${PROJECT_LIBS_FLAGS} ${MPI_LIB} ${PROJECT_LIB_DIR_FLAG} ${STATIC}" \
   MPI=$(DEF_V_MPI) DEFINES="${DEFINE4ALL} -D_DEBUG" $(PN)
#========================================================================================================#

#=================================== vlpl3d project SINGLE VERSION ======================================#
single : _single strip

_single :
	$(M_MAKE) OO=$(Optimized) LF="${PROJECT_LIBS_FLAGS} ${PROJECT_LIB_DIR_FLAG} ${STATIC}" DEFINES="${DEFINE4ALL}" \
   $(PN)
#========================================================================================================#

#================================ vlpl3d project SINGLE VERSION DEBUG====================================#
singled : _singled

_singled :
	$(M_MAKE) OO=$(Debug) GG=-g LF="${PROJECT_LIBS_FLAGS} ${PROJECT_LIB_DIR_FLAG} ${STATIC}" \
   DEFINES="${DEFINE4ALL} -D_DEBUG" $(PN)
#========================================================================================================#

#=========================== remove debug information from release file =================================#           
strip :
	strip -s $(PN)
   
strippp :
	strip -s $(PPN)
#========================================================================================================#
   
#===================================== vlpl3d project CLEAN =============================================#           
allclean : exec_clean clean libsclean

exec_clean :
	-rm -f *.e

clear : clean

clean :     
	-rm -f *.o

libsclean : pulse_profile_libs_clean scales_libs_clean projects_libs_clean   
      
pulse_profile_libs_clean :
	$(MAKE) -C ${PULSE_PROFILE_LIB_DIR} clean   

scales_libs_clean :
	$(MAKE) -C ${SCALES_LIB_DIR} clean
      
projects_libs_clean :
	-rm -f $(PROJECT_LIB_DIR)/*.a

.PHONY : clean
#========================================================================================================#

#======================================= objects compile ================================================#              
.C.o	  :
		$(VLPLCC) $(CFLAGS) $<

.cpp.o	  :
		$(VLPLCC) $(CFLAGS) $<

.c.o	  :
		$(VLPLCC) -Kc++ $(CFLAGS) $<
#========================================================================================================#      

#======================================= Compile on Mac OS X (rsg) ======================================#

macosx :
	mpic++ -I/opt/local/include -L/opt/local/lib *.cpp -lhdf5 -lz -o vdsr_macosx.e

