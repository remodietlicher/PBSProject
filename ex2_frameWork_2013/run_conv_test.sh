#!/bin/bash

file=log.log;

echo "Convergence Test: ";
for i in 3 5 8 16 32;
do
    echo ${i};
    echo GRISIZE = $i >> $file;
    clang++ *.cpp -framework OpenGL -lglut -DGRIDSIZE=${i} && ./a.out | grep Error >> $file;
done
echo "results in: $file"
