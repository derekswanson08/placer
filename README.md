# Placer

A quadratic placer for laying out gates on an integrated circuit, or components on a printed circuit board.

The input is a netlist describing how the gates connect to each other, as well as how the gates connect to pads at fixed locations around the perimeter of the chip. The output is the X/Y position of each gate on the chip surface.

The optimal solution is defined as the placement of all gates such that the total wirelength is minimized. This typically results in a lot of gates being placed near or on top of eachother. Partitioning is used to more evenly spread the gates out across the surface of the chip.

At each partitioning step, a subregion of the chip is divided into two halves, and half the gates are assigned to each half. Then the quadratic placement is applied locally to each half. This procedure is iterated until the gates are sufficiently spread out.

## Build/Run
* Open placer.pro in Qt Creator
* Click Build/Run

## Dependencies
* [Qt 5.14](https://www.qt.io) GUI library
* [Eigen](http://eigen.tuxfamily.org) Linear algebra library

