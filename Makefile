CFLAGS = -Drestrict=__restrict__ -DREUSE_CSR_FOR_VALIDATION -I../aml \
	-march=znver4 -Ofast -flto -funroll-loops -fomit-frame-pointer -fno-math-errno -DNDEBUG -fopenmp

# Optionnel: pointer TBB_DIR (ex: export TBB_DIR=$(spack location -i tbb))
TBB_DIR ?=
TBB_INC = $(if $(TBB_DIR),-I$(TBB_DIR)/include,)
TBB_LIB = $(if $(TBB_DIR),-L$(TBB_DIR)/lib -L$(TBB_DIR)/lib64,)

CFLAGS += $(TBB_INC)
LDFLAGS = -lm -lpthread -flto $(TBB_LIB) -ltbb -ltbbmalloc -ltbbmalloc_proxy
MPICC = mpicxx

all: main

GENERATOR_SOURCES = generator/graph_generator.cpp generator/make_graph.cpp generator/splittable_mrg.cpp generator/utils.cpp generator/generator.cpp kernels/gen_graph.cpp kernels/shortest_path.cpp kernels/breadth_search.cpp
SOURCES = main.cpp
HEADERS = 

main: $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CFLAGS) -o main $(SOURCES) $(GENERATOR_SOURCES) $(LDFLAGS) $(CXXFLAGS) $(LDLIBS)