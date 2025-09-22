# KFParticle.hxx

A fork of [KFParticle](https://github.com/alisw/KFParticle), with:

- same mathematics
- more protections
- slightly faster, with cached output of expensive operations
- reduced to contain only ALICE-related settings (homogeneous magnetic field in z-axis)
- formatted with `.clang-format` and `.clang-tidy`
- more comments
- upgraded CMake instructions
- no third-party library dependencies

## Requirements

- CMake (v3.25 or higher)
- C++ compiler compatible with C++23

## Build

```bash
mkdir <build-dir> && cd <build-dir>
cmake <source-dir> <options>
cmake --build .
cmake --install . --prefix <install-dir>
```

Additional options:

* `-DKF_DEBUG=ON` -- enable debug messages

## Test

Using `g++`:

```bash
cd KFParticleTest/
g++ KFParticleTest.cxx -std=c++23 -march=native -mtune=native -O3 -DNDEBUG -I<install-dir>/include -L<install-dir>/lib -lKFParticle -o Test
./Test # output should be equivalent to KFParticleTest/test.txt
```

## CMake Integration

Add these lines or equivalent to your `CMakeLists.txt` files:

```cmake
message(STATUS "Looking for KFParticle")
find_package(KFParticle)
if(KFParticle_FOUND)
    message(STATUS "Looking for KFParticle -- found")
    message(STATUS "KFParticle found at ${KFParticle_LIB_DIR}")
endif()

target_link_libraries(App PUBLIC KFParticle::KFParticle)
```

And make sure to prepare your CMake cache (step before building) with `-DKFParticle_DIR=<kfparticle-install-dir>/lib/cmake/KFParticle/`.

## Examples

* [tree2secondaries](https://github.com/HD-ALICE-Sexaquark/tree2secondaries)
