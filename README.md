# SA-f(PSI)

install required packages: sudo apt-get install libboost-all-dev

How to compile and run:

1. g++ -O3 -ftree-vectorize -march=native -msse4.1 -o PSI PSI.cpp Circuit/Gate.cpp Utility/ISecureRNG.cpp Utility/CryptoUtility.cpp Utility/Commitment.cpp Utility/Range.cpp Utility/Communicator.cpp Polynomial.cpp -lssl -lntl -lgmp -lpthread -std=c++11 -lcrypto -l:/usr/local/lib/libcryptopp.a -maes -pthread -mrdseed -rdynamic /usr/local/lib/librelic.so -lboost_system -lgmp /usr/local/lib/libemp-tool.so -Wl,-rpath,/usr/local/lib

Link emp example:

g++ -std=c++11 shot.cpp -pthread -Wall -march=native -O3 -maes -mrdseed -rdynamic /usr/local/lib/librelic.so -lssl -lcrypto -lboost_system -lgmp /usr/local/lib/libemp-tool.so -Wl,-rpath,/usr/local/lib

2. $ ./PSI <number of parallel machines> <machine Id> <parallel peers base port> <party Id> <parties base port> <input length> <test no> <|PSI| fraction>

You can modify machine ip configurations in machine_spec folder.

clear && ./PSI 1 0 40000 0 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_1.txt"
clear && ./PSI 1 0 40000 0 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_2.txt"
clear && ./PSI 1 0 40000 0 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_3.txt"
clear && ./PSI 1 0 40000 0 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_4.txt"
clear && ./PSI 1 0 40000 0 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_5.txt"

clear && ./PSI 1 0 40000 1 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_1.txt"
clear && ./PSI 1 0 40000 1 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_2.txt"
clear && ./PSI 1 0 40000 1 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_3.txt"
clear && ./PSI 1 0 40000 1 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_4.txt"
clear && ./PSI 1 0 40000 1 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_5.txt"

clear && ./PSI 1 0 40000 2 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_1.txt"
clear && ./PSI 1 0 40000 2 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_2.txt"
clear && ./PSI 1 0 40000 2 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_3.txt"
clear && ./PSI 1 0 40000 2 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_4.txt"
clear && ./PSI 1 0 40000 2 50000 1048576 1 0.25 |& tee "nn20_circuit_psi_0.25_5.txt"

clear && ./PSI 1 0 40000 3 50000 1048576 3 0.5

./FourParty $1 $i 40000 $2 50000 $3 $4 $5 0.3 |& tee "$i.txt" & #|& tee "q$i.txt" &
done
