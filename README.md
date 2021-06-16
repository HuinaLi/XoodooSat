
# Xoodoo SAT solver

**Using [CryptominiSAT](https://github.com/msoos/cryptominisat/) to find differential trails of [Xoodoo](https://keccak.team/xoodoo.html)**

## Install

**1. Install cryptominisat**

To build and install, first get the tar.gz package in [here](https://github.com/msoos/cryptominisat/releases), and issue(in Linux):


```
sudo apt-get install build-essential cmake
# not required but very useful
sudo apt-get install zlib1g-dev libboost-program-options-dev libm4ri-dev libsqlite3-dev help2man
tar xzvf cryptominisat-version.tar.gz
cd cryptominisat-version
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```

**2. Install pysat**

We use [pysat](https://github.com/pysathq/pysat) to generate the **cnf file** for cardinality constraint, which is used to bound the weight of the trail.
A cardinality constraint is a constraint of the form: $$\sum_{i=1}^n{x_i}\leq k, x_i={0,1}$$ Here n is the number of variables.

```
pip3 install python-sat
```

**3. Compile XOODOOSAT**

Just clone and run:

```
make
```

By default the output executable file is **xoodoo** defined in **Makefile**.

## Test

Run:

```
./xoodoo -h
```
you can see the optional parameters like round number(how many rounds to analysis), weight, etc.

An example(3 rounds, max weight 25):
```
./xoodoo -r 3 -w 25 -t 16 -m 0
```

Note that:

$$3 round: a_0 -> b_0 -> a_1 -> b_1 -> a_2 -> b_2 -> a_3, b_i=lambda(a_i), a_i+1=chi(b_i)$$
$$3 round trail core: a_1 -> b_1 -> a_2 -> b_2$$
$$state to sum weight: a_1, a_2, b_2 (a_1, a_2,..., a_n-1, b_n-1)$$
$$so 3 round weight = w(a_1)+w(a_2)+w(b_2)$$


Finally, the result is output in result folder.