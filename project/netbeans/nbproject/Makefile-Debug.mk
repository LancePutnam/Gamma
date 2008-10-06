#
# Gererated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Ambisonics.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../tests/testAudioChannels.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Delay.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/DataFile.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/arr.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/scl.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Sync.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/SoundFile.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/AudioIO.o \
	${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Oscillator.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-lpthread -ldl /usr/lib/libportaudio.a /usr/lib/libsndfile.a /usr/lib/libjack.a /usr/lib/libFLAC.a /usr/lib/libogg.a /usr/lib/libasound.a /usr/lib/libc_nonshared.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} dist/Debug/GNU-Linux-x86/netbeans

dist/Debug/GNU-Linux-x86/netbeans: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -lrt -o dist/Debug/GNU-Linux-x86/netbeans ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Ambisonics.o: ../src/Ambisonics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Ambisonics.o ../src/Ambisonics.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../tests/testAudioChannels.o: ../tests/testAudioChannels.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../tests
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../tests/testAudioChannels.o ../tests/testAudioChannels.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Delay.o: ../src/Delay.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Delay.o ../src/Delay.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/DataFile.o: ../src/DataFile.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/DataFile.o ../src/DataFile.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/arr.o: ../src/arr.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/arr.o ../src/arr.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/scl.o: ../src/scl.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/scl.o ../src/scl.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Sync.o: ../src/Sync.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Sync.o ../src/Sync.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/SoundFile.o: ../src/SoundFile.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/SoundFile.o ../src/SoundFile.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/AudioIO.o: ../src/AudioIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/AudioIO.o ../src/AudioIO.cpp

${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Oscillator.o: ../src/Oscillator.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src
	$(COMPILE.cc) -g -I../include -o ${OBJECTDIR}/_ext/home/lance/projects/libs/synz/netbeans/../src/Oscillator.o ../src/Oscillator.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/netbeans

# Subprojects
.clean-subprojects:
