CFLAGS = -Drestrict=__restrict__ -DREUSE_CSR_FOR_VALIDATION -I../aml \
	-march=znver4 -Ofast -flto -funroll-loops -fomit-frame-pointer -fno-math-errno -DNDEBUG -fopenmp
LDFLAGS = -lm -lpthread -ltbb -flto
MPICC = mpicxx

all: main

GENERATOR_SOURCES = generator/graph_generator.cpp generator/make_graph.cpp generator/splittable_mrg.cpp generator/utils.cpp generator/generator.cpp kernels/gen_graph.cpp kernels/shortest_path.cpp kernels/breadth_search.cpp
SOURCES = main.cpp
HEADERS = 

main: $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CFLAGS) -o main $(SOURCES) $(GENERATOR_SOURCES) $(LDFLAGS) $(CXXFLAGS) $(LDLIBS)