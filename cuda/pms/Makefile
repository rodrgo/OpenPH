include ../config  

#VERSION=$(shell date +"%Y%m%d")
VERSION="_0_0_0"


TARGET_SMV=ph.mexa64 

TARGET_PACKAGE=TDA$(VERSION).tar

.PHONY: all smv

all: $(TARGET_SMV)

clean:
	rm -f $(TARGET_SMV)

smv: $(TARGET_SMV)

package: $(TARGET_PACKAGE)

clean_package:
	rm -f ../../$(TARGET_PACKAGE)

%.mexa64: %.cu
	$(MAKE) -f $(CUDAMATLAB)/Makefile $@

