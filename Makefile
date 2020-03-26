FLAGS = -m64 -std=c++17 -Wall -Werror -Wextra -pedantic -Wformat -Wformat-security -Wstrict-overflow -Wconversion -ftree-vectorize -O2 -g
SIMPLEX:
	g++ *.hpp *.cpp -o main $(FLAGS)