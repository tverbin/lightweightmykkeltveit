## Mykkeltveit's Decycling Algorithm Implementation
This is an efficient lightweight C implementation to the Mykkeltveit's Algorithm for decycling de-Bruijn graph.
The alphabet used for the code is '0123', but can be replaced to any 2^k long alphabet.

#### Usage

    make
    ./decycling <k>

#### Credits
The concept is based on Mykkeltveit's article from 1971:
Mykkeltveit J. A proof of Golomb’s conjecture for the de Bruijn graph. Journal of Combinatorial Theory, Series B. 1972 Aug;13(1):40–45.
The implementation is inspired by a java implementation done in Ron Shamir Laboratory in Tel Aviv University:
https://github.com/shamir-lab/docks

#### Licence
GPLv2
