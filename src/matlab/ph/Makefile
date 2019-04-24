include ../config  

#VERSION=$(shell date +"%Y%m%d")
VERSION="_0_0_0"

TARGET_SMV=ph.mexa64 

.PHONY: all smv

all: $(TARGET_SMV)

clean:
	rm -f $(TARGET_SMV)

smv: $(TARGET_SMV)

%.mexa64: %.cu
	$(MAKE) -f $(CUDAMATLAB)/Makefile $@

