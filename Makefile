OBJS = make_coefficient_matrix.o  initialize.o set_array.o inner_product.o set_variables.o  gaussian_pulse.o main.o pml.o tfsf.o dielectrics.o ununiform.o interpolation.o fourie_transform.o
OPTS = -O3 -std=c++1y -I/opt/include/eigen3/
HEAD = make_coefficient_matrix.h  initialize.h set_array.h inner_product.h set_variables.h gaussian_pulse.h pml.h tfsf.h dielectrics.h ununiform.h interpolation.h fourie_transform.h

all:CIPK 

.PHONY: all 

CIPK: $(OBJS) $(HEAD)
	g++ -o $@.exe $(OBJS) $(OPTS) 

%.o: %.cpp $(HEAD)
	g++ -c $< $(OPTS)

