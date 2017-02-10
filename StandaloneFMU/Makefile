# Makefile for 20-sim FMU %FMIVERSION% export

FMU_NAME=%SUBMODEL_NAME%
FMU=$(FMU_NAME).fmu

FMU_DIR=project/fmu
FMU_SRC_DIR   = $(FMU_DIR)/sources
FMU_DOC_DIR   = $(FMU_DIR)/documentation

BUILD = project/gcc
SRC   = src

XXSIM_SOURCES = EulerAngles.c \
	MotionProfiles.c \
	xxfuncs.c \
	xxinteg.c \
	xxinverse.c \
	xxmatrix.c \
	xxmodel.c \
	xxsubmod.c \
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
	cvode/src/cvode/cvode.c \
	cvode/src/cvode/cvode_dense.c \
	cvode/src/cvode/cvode_direct.c \
	cvode/src/cvode/cvode_io.c \
	cvode/src/nvec_ser/nvector_serial.c \
	cvode/src/sundials/sundials_dense.c \
	cvode/src/sundials/sundials_direct.c \
	cvode/src/sundials/sundials_math.c \
	cvode/src/sundials/sundials_nvector.c \
%ENDIF%
%IF%%EQ(INTEGRATION_METHOD_NAME,MeBDFi)%
	MeBDFi/MeBDFi.c \
	MeBDFi/MeBDFi.h \
	MeBDFi/include/f2c.h \
	MeBDFi/include/fio.h \
	MeBDFi/include/fmt.h \
	MeBDFi/include/fp.h \
	MeBDFi/include/lio.h \
	MeBDFi/include/sysdep1.h \
	MeBDFi/libf2c/close.c \
	MeBDFi/libf2c/dolio.c \
	MeBDFi/libf2c/endfile.c \
	MeBDFi/libf2c/err.c \
	MeBDFi/libf2c/fmt.c \
	MeBDFi/libf2c/fmtlib.c \
	MeBDFi/libf2c/lread.c \
	MeBDFi/libf2c/lwrite.c \
	MeBDFi/libf2c/open.c \
	MeBDFi/libf2c/pow_dd.c \
	MeBDFi/libf2c/pow_di.c \
	MeBDFi/libf2c/rdfmt.c \
	MeBDFi/libf2c/rsfe.c \
	MeBDFi/libf2c/sfe.c \
	MeBDFi/libf2c/sig_die.c \
	MeBDFi/libf2c/util.c \
	MeBDFi/libf2c/wref.c \
	MeBDFi/libf2c/wrtfmt.c \
	MeBDFi/libf2c/wsfe.c \
	MeBDFi/libf2c/wsle.c \
	MeBDFi/netlib/DOGLEG.c \
	MeBDFi/netlib/DPMPAR.c \
	MeBDFi/netlib/ENORM.c \
	MeBDFi/netlib/FDJAC1.c \
	MeBDFi/netlib/HYBRD.c \
	MeBDFi/netlib/HYBRD1.c \
	MeBDFi/netlib/mebdfiCore.c \
	MeBDFi/netlib/QFORM.c \
	MeBDFi/netlib/QRFAC.c \
	MeBDFi/netlib/R1MPYQ.c \
	MeBDFi/netlib/R1UPDT.c \
%ENDIF%
	xxTable2D.c

FMU_SOURCES = %FMI_PREFIX%Functions.c $(XXSIM_SOURCES)

INCLUDES=
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
INCLUDES += -Isrc/cvode/include

%ENDIF%
%IF%%EQ(INTEGRATION_METHOD_NAME,MeBDFi)%
INCLUDES += -Isrc/MeBDFi/include

%ENDIF%
OBJS = $(addprefix $(BUILD)/, $(addsuffix .o, $(basename $(FMU_SOURCES))))
DEPS = $(addprefix $(BUILD)/, $(addsuffix .d, $(basename $(FMU_SOURCES))))

CXX?=g++
CC?=gcc

current_dir = $(shell pwd)

# OS detection
ifeq ($(OS),Windows_NT)
	# Windows
	WINDOWS=1
	DEFINES=-DWIN32
	# For now, assume 32-bit
	# TODO detect MinGW 64
	OS_BITS=32
	DLL_EXT=dll
	BUILD_DLL=-shared
else
	# Other OS (Linux/MacOSX)
	ifeq ($(shell uname), Linux)
		LINUX=1
		DEFINES=-DLINUX
		DLL_EXT=so
		BUILD_DLL=-shared
	endif
	ifeq ($(shell uname), Darwin)
		MACOSX=1
		DEFINES=-DMACOSX
		DLL_EXT=so
		BUILD_DLL=-dynamiclib
	endif
	# Detect 32-bit / 64-bit OS
	OS_BITS=$(shell getconf LONG_BIT)
endif

# Set the OS specific settings
ifeq ($(WINDOWS),1)
	FMU_BIN32_DIR = $(FMU_DIR)/binaries/win32
	FMU_BIN64_DIR = $(FMU_DIR)/binaries/win64
	TEMPLATE_DIR=$(subst :,,$(subst \,/,/%TEMPLATE_DIR%))
	ZIP_COMMAND=$(TEMPLATE_DIR)/bin/7z.exe a -tzip -xr!.svn "$(FMU)" *
	POS_INDEP_CODE=
endif
ifeq ($(LINUX),1)
	FMU_BIN32_DIR = $(FMU_DIR)/binaries/linux32
	FMU_BIN64_DIR = $(FMU_DIR)/binaries/linux64
	ZIP_COMMAND=zip -r $(FMU) *
	POS_INDEP_CODE=-fPIC
endif
ifeq ($(MACOSX),1)
	FMU_BIN32_DIR = $(FMU_DIR)/binaries/darwin32
	FMU_BIN64_DIR = $(FMU_DIR)/binaries/darwin64
	ZIP_COMMAND=zip -r $(FMU) *
	POS_INDEP_CODE=-fPIC
endif

ifeq ($(MAKECMDGOALS),debug)
  CFLAGS=-g -DDEBUG $(INCLUDES) $(DEFINES) $(POS_INDEP_CODE)
  DEBUG=1
else
  CFLAGS=-O2 -Wall $(INCLUDES) $(DEFINES) $(POS_INDEP_CODE)
  DEBUG=0
endif

debug: fmu

default: fmu

.PHONY: clean

clean:
	@rm -rf $(BUILD)
	@rm -f $(FMU)
	@rm -rf $(FMU_BIN32_DIR)
	@rm -rf $(FMU_BIN64_DIR)
	@rm -rf $(FMU_SRC_DIR)
	@rm -f $(FMU_DOC_DIR)/*

welcome:
	@echo ------------------------------------------------------------
	@echo 20-sim standalone co-simulation FMU export for:
	@echo '$(FMU_NAME)'
	@echo ------------------------------------------------------------

$(BUILD):
	@mkdir -p $(BUILD)
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
	@mkdir -p $(BUILD)/cvode/src/cvode
	@mkdir -p $(BUILD)/cvode/src/nvec_ser
	@mkdir -p $(BUILD)/cvode/src/sundials
%ENDIF%
%IF%%EQ(INTEGRATION_METHOD_NAME,MeBDFi)%
	@mkdir -p $(BUILD)/MeBDFi/include
	@mkdir -p $(BUILD)/MeBDFi/libf2c
	@mkdir -p $(BUILD)/MeBDFi/netlib
%ENDIF%

$(BUILD)/%.o: $(SRC)/%.c Makefile $(BUILD)
	@echo "[$(CC)] $<"
	@$(CC) -MD -MF $(patsubst %.o,%.d,$@) -c $(CFLAGS) $< -o $@

dll: welcome $(OBJS)
	@$(CC) $(CFLAGS) -o $(BUILD)/$(FMU_NAME).$(DLL_EXT) $(BUILD_DLL) $(OBJS) -lm
	@[ $(DEBUG)==0 ] && strip $(BUILD)/$(FMU_NAME).$(DLL_EXT)

createfmudir:
	@echo Creating an empty FMU
	@mkdir -p $(FMU_DIR)
	@mkdir -p $(FMU_SRC_DIR)
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
	@mkdir -p $(FMU_SRC_DIR)/cvode
	@mkdir -p $(FMU_SRC_DIR)/nvector
	@mkdir -p $(FMU_SRC_DIR)/sundials
%ENDIF%
fmu: dll createfmudir
	@if [ $(OS_BITS) -eq 32 ]; then \
		mkdir -p $(FMU_BIN32_DIR); \
		echo Copy the compiled 32-bit $(DLL_EXT) $(BUILD)/$(FMU_NAME).$(DLL_EXT) to $(FMU_BIN32_DIR); \
		cp $(BUILD)/$(FMU_NAME).$(DLL_EXT) $(FMU_BIN32_DIR); \
	else \
		mkdir -p $(FMU_BIN64_DIR); \
		echo Copy the compiled 64-bit $(DLL_EXT) $(BUILD)/$(FMU_NAME).$(DLL_EXT) to $(FMU_BIN64_DIR); \
		cp $(BUILD)/$(FMU_NAME).$(DLL_EXT) $(FMU_BIN64_DIR); \
	fi

	@echo Copy the generated sources to $(FMU_SRC_DIR)
	@cp $(SRC)/*.c $(FMU_SRC_DIR)
	@cp $(SRC)/*.h $(FMU_SRC_DIR)
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
	@cp $(SRC)/includes.txt $(FMU_SRC_DIR)
	@cp $(SRC)/cvode/src/cvode/cvode.c $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/cvode/cvode_dense.c $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/cvode/cvode_direct.c $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/cvode/cvode_io.c $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/nvec_ser/nvector_serial.c $(FMU_SRC_DIR)/nvector/
	@cp $(SRC)/cvode/src/sundials/sundials_dense.c $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/src/sundials/sundials_direct.c $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/src/sundials/sundials_math.c $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/src/sundials/sundials_nvector.c $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/cvode/cvode.h $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/include/cvode/cvode_direct.h $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/include/cvode/cvode_dense.h $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/cvode/cvode_direct_impl.h $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/src/cvode/cvode_impl.h $(FMU_SRC_DIR)/cvode/
	@cp $(SRC)/cvode/include/nvector/nvector_serial.h $(FMU_SRC_DIR)/nvector/
	@cp $(SRC)/cvode/include/sundials/sundials_math.h $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/sundials/sundials_nvector.h $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/sundials/sundials_config.h $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/sundials/sundials_dense.h $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/sundials/sundials_direct.h $(FMU_SRC_DIR)/sundials/
	@cp $(SRC)/cvode/include/sundials/sundials_types.h $(FMU_SRC_DIR)/sundials/
%ENDIF%
%IF%%EQ(INTEGRATION_METHOD_NAME,MeBDFi)%
	@cp $(SRC)/includes.txt $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/MeBDFi.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/MeBDFi.h $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/include/f2c.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/include/fio.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/include/fmt.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/include/fp.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/include/lio.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/include/sysdep1.h $(FMU_SRC_DIR)/include
	@cp $(SRC)/MeBDFi/libf2c/close.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/dolio.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/endfile.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/err.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/fmt.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/fmtlib.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/lread.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/lwrite.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/open.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/pow_dd.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/pow_di.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/rdfmt.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/rsfe.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/sfe.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/sig_die.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/util.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/wref.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/wrtfmt.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/wsfe.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/libf2c/wsle.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/DOGLEG.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/DPMPAR.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/ENORM.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/FDJAC1.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/HYBRD.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/HYBRD1.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/mebdfiCore.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/QFORM.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/QRFAC.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/R1MPYQ.c $(FMU_SRC_DIR)
	@cp $(SRC)/MeBDFi/netlib/R1UPDT.c $(FMU_SRC_DIR)
%ENDIF%
	@echo Create the FMU: $(FMU)
	@cd $(FMU_DIR) && $(ZIP_COMMAND) && mv $(FMU) $(current_dir) && echo Ready!

