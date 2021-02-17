all:
	clang algorithm.c pdb.c -O0 -g -fsanitize=address -pedantic-errors -lm -o align
	chmod +x check_3dna.sh
