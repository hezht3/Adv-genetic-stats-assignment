cd /home/Johnathan/caviar/CAVIAR-C++
./CAVIAR -h
./CAVIAR -c 1 -l input/LD.txt -z input/Z.txt -o out/out_causal_1 -f 1
./CAVIAR -c 2 -l input/LD.txt -z input/Z.txt -o out/out_causal_2 -f 1
./CAVIAR -c 3 -l input/LD.txt -z input/Z.txt -o out/out_causal_3 -f 1

./CAVIAR -c 1 -l input/LD.txt -z input/Z.txt -r 0.95 -o out/out_causal_1_95 -f 1
./CAVIAR -c 2 -l input/LD.txt -z input/Z.txt -r 0.95 -o out/out_causal_2_95 -f 1
./CAVIAR -c 3 -l input/LD.txt -z input/Z.txt -r 0.95 -o out/out_causal_3_95 -f 1

./CAVIAR -c 1 -l input/LD.txt -z input/Z.txt -r 0.99 -o out/out_causal_1_99 -f 1
./CAVIAR -c 2 -l input/LD.txt -z input/Z.txt -r 0.99 -o out/out_causal_2_99 -f 1
./CAVIAR -c 3 -l input/LD.txt -z input/Z.txt -r 0.99 -o out/out_causal_3_99 -f 1

