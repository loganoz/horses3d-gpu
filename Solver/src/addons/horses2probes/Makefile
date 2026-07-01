#####################################################################
## HORSES - horses2probes addon Makefile
#####################################################################

INSTALL_DIR ?= ./
PROG        = ./bin/horses2probes

# IMPORTANT: force FLOW for this addon
# Use 'override' because the top-level Makefile invokes this sub-make with
# MACROS='$(MACROS)' on the command line, which would otherwise take
# precedence over a plain '+=' assignment here and silently drop -DFLOW.
override MACROS += -DFLOW

OBJS =  sharedSpectralBasis \
        Storage \
        Solution2probes \
        main

OBJS_BUILD = $(foreach obj,$(OBJS),./build/$(obj).o)

-include ./make.deps

.DEFAULT_GOAL := $(PROG)

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) $(INCLUDE)_ns \
	      -c $< -o $@ $(MOD_FLAG) ./include

$(PROG): header mkdirs $(OBJS_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(EXTLIB_INC) -o $(PROG) $(OBJS_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)

install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
	cp -p ./include/*.mod $(INSTALL_DIR)/include
	@echo
	@echo

clean: FORCE
	rm -f ./include/*.mod ./build/*.o
	rm -f $(PROG)

mkdirs: FORCE
	@mkdir -p ./build
	@mkdir -p ./include
	@mkdir -p ./bin

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building horses2probes"
	@echo "================================"
	@echo
	@echo

FORCE:
