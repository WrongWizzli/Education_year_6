g++ solver.cpp -O3 -fopenmp -osolver -std=c++11

nthreads=(1 2 4 8 16 32 64 128)
L=(1 3.14159265358979323846264338327950288)
space_grid=(64 128 256 512)
T=(0.05 0.025)
K=(100 40)

for i in ${!nthreads[@]}; do
    for j in ${!L[@]}; do
        for k in ${!space_grid[@]}; do
            echo "bsub ./solver ${L[$j]} ${L[$j]} ${L[$j]} 0.05 ${space_grid[$k]} 100 ${nthreads[$i]}"
            echo "bsub ./solver ${L[$j]} ${L[$j]} ${L[$j]} 0.05 ${space_grid[$k]} 100 ${nthreads[$i]}" >> launches.txt
            bsub ./solver ${L[$j]} ${L[$j]} ${L[$j]} 0.05 ${space_grid[$k]} 100 ${nthreads[$i]}
        done
    done
done