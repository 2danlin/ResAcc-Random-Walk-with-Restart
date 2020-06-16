g++ -march=core2 -ffast-math -std=c++11 -pthread -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o ResAcc SFMT.c main.cpp
