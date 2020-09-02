# Constant-time ROLLO

A constant time implementation of ROLLO_I_128 with AVX2 instructions.

This project contains the code described in the paper:

_Constant time algorithms for ROLLO-I-128_

## Authors
- Carlos Aguilar-Melchor, 
  ISAE-SUPAERO, Université de Toulouse, Toulouse, France 
- Nicolas Aragon, Université de Limoges, Limoges Cedex, France
- Emanuele Bellini, Cryptography Research Centre, Technology Innovation Institute, Abu Dhabi, UAE 
- Florian Caullery , Cryptography Research Centre, Technology Innovation Institute, Abu Dhabi, UAE
- Rusydi H. Makarim, Cryptography Research Centre, Technology Innovation Institute, Abu Dhabi, UAE
- Chiara Marcolla, Cryptography Research Centre, Technology Innovation Institute, Abu Dhabi, UAE

## Compilation instructions

The project can either be run using CLion or from the command line.

### To run from command line

From a terminal, follow the steps:

```
git clone https://github.com/peacker/constant_time_rollo.git
cd constant_time_rollo
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release -DOPTIMIZATION=AVX2
make
./rollo-i-128/rollo_i_128_simple_test
./tests/rollo_i_128_tests
```

The command

`./rollo-i-128/rollo_i_128_simple_test`

simply runs Key generation, Encapsulation and Decapsulation.

The command

`./tests/rollo_i_128_tests`

runs unit and performance tests.
