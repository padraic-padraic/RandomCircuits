platform = $(shell uname)
pybind_includes = $(shell python3 -m pybind11 --include)
include_dirs = -I./ -I./circuit_generation $(pybind_includes)
pybind_extension = $(shell python3-config --extension-suffix)
CCX = g++
CPP_FLAGS = -Wall -O3 -shared -std=c++14
pybind_flags = -fPIC
pybind_darwin = -undefined dynamic_lookup

nlohmann = ./nlohmann/json.hpp

all: extension

download:
ifeq ("$(wildcard $(nlohmann))","")
	@echo Downloaded nlohmann json
	$(shell curl https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp -o nlohmann/json.hpp)
else
	@echo Found nlohmann json, skipping the download step...
endif

ifeq ($(platform), Darwin)
CPP_FLAGS += $(pybind_darwin)
else
CPP_FLAGS += $(pybind_flags)
endif

source = ./circuit_generation/random_circuits.cpp
target := $(addsuffix ${pybind_extension}, ./circuit_generation/random_circuits)

extension: download $(source)
	$(CCX) $(CPP_FLAGS) $(include_dirs) $(source) -o $(target)

.PHONY: clean
clean:
	rm $(target)

