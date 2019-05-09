
INC_DIR := -I/usr/local/opt/petsc/include
LIB_DIR := -L/usr/local/opt/petsc/lib
#LIB_DIR := -L/home/alain/Bureau/Projet_2660/lib_petsc/lib -Wl,-rpath=/home/alain/Bureau/Projet_2660/lib_petsc/lib
LIB := -lpetsc.3.10

CXX_FLAGS := -O0 -Wall #-g

#Compilation 
all : 
	gcc -o projet projet.c Write.c ns_solver.c poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Delete of executable files
clean :
	rm projet

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
