# TFM-Quantum-Kernel-SVM

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)


## General info
This project's motivation is to study the viability of quantum computing in machine learning systems's training proccess. We'll study its use in classification problems and Suppor Vector Machines. 

We'll use the work developed by Havlíček, V., Córcoles, A.D., Temme, K. et al in their publication ["Supervised learning with quantum-enhanced feature spaces". Nature 567, 209–212 (2019](https://arxiv.org/pdf/1804.11326.pdf) as the baseline for our studies.

In our work we'll focus on how dimensionality can affect the time required to train a SVM classifier. For this purpouse we'll train an SVM using both classical and quantum kernels an compare the change in training times when we increase dimensionality.

We spect a more robust behavior by the quantum kernel when we have high dimensionality. 

## Technologies
Project is created with:
* Python version: 3.X
* Qiskit version: XX.XX
* Qiskit_Machine_Learning: XX.XX
