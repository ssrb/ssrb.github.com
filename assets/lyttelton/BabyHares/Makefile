V = 0

AT_0 := @
AT_1 :=
AT = $(AT_$(V))

REDIR_0 := > /dev/null
REDIR_1 :=
REDIR = $(REDIR_$(V))

CXX := /opt/intel/bin/icpc
CXXFLAGS := -O3 -I/usr/lib/openmpi/include \
	-I/usr/lib/openmpi/include/openmpi \
	-I/usr/local/cuda/include \
	-mkl=parallel
# -llmpe -lmpe
LDFLAGS :=  -pthread -openmp -mkl=parallel \
	-L/usr/lib/openmpi/lib -lmpi -lopen-rte -lopen-pal\
	-L/usr/local/cuda/lib64 \
	-L/usr/local/cuda/lib \
	-ldl -lnsl -lutil \
	-lcudart \
	-lcusparse 

NVCC := /usr/local/cuda/bin/nvcc
NVCCFLAGS := -O3 -arch=sm_20 -Xcompiler -Wall -Xcompiler -Wextra -m64

CPPSRC=$(shell find . -name "*.cpp")
CUDASRC=$(shell find . -name "*.cu")

OUT=obj
OBJS=$(CPPSRC:%.cpp=$(OUT)/%.o) $(CUDASRC:%.cu=$(OUT)/%.o) 
DEPS=$(CPPSRC:%.cpp=$(OUT)/%.d)

TARGET=babyharesED

NODEPS=clean

.PHONY = all clean

.SECONDEXPANSION:

all: $(TARGET)

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
    -include $(DEPS)
endif

$(TARGET): $(OBJS)
	@echo [LD] $@
	$(AT)$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@ $(REDIR)

$(OUT)/%.d:%.cpp $$(@D)/.f
	@echo [DEP] $<
	$(AT)$(CXX) $(CXXFLAGS) -MM -MT '$(patsubst %.d,%.o,$@)' $< -MF $@ $(REDIR)

$(OUT)/%.o:%.cpp $(OUT)/%.d
	@echo [C++] $<
	$(AT)$(CXX) $< $(CXXFLAGS) -c -o $@ $(REDIR)

$(OUT)/%.o:%.cu
	@echo [CUDA] $<
	$(AT)$(NVCC) $(NVCCFLAGS) $< -c -o $@ $(REDIR)

%/.f:
	@echo [MKDIR] $(dir $@)
	$(AT)mkdir -p $(dir $@) 
	$(AT)touch $@

.PRECIOUS: %/.f

clean:
	@echo [CLEAN]
	$(AT)-rm -rf $(TARGET) $(OUT)
