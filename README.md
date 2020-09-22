# Fox-Algorithm-Parallel-Matrix-Multiply


>Implementation of Fox' algorithm for matrix multiplication.
Based on Matrix-Matrix Multiplication Using Fox’s Algorithm and MPI standard.

## Table of contents
* [About Fox Algorithm](#about-fox-algorithm)
* [Technologies](#technologies)
* [Status](#status)
* [Contact](#contact)

## About Fox Algorithm
We are going to have n x n tasks
Each task will manage the appropriate block entry in the matrices.
The entries in the matrix are matrices of size blk.
Thus in reality we will be multiplying n x blk square matrices.

>How the algorithm works?

Each task starts out with an A , B , and C block.
In step one, the tasks on the diagonal  broadcast their  block along their row.
They each then multiply this A block with their B block and add it to their C block.
Now we rotate the B blocks.
We pass our B block to the task above us (if we are in the top row we pass our B to the entry in the bottom row) but if you think of the thing as a torus we are just passing it up.


## Technologies

* C
* MPI

## Status
Project is: _finished_

## Contact
Created by [@ValentinaSimic](https://github.com/ValentinaSimic) - feel free to contact me!
