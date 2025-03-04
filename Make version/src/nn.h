#ifndef NN_H
#define NN_H
#include <random>
#include <cmath>
#include "matrix.h"

extern float h; //step to calculate derivatives
extern float u; //regularization parameter
extern float u_max; //max value of 'u'
extern float u_min; //min value of 'u'
extern float u_inc; //coefficient to increase 'u'
extern float u_dec; //coefficient to decrease 'u'
extern float err_min; //error when stop
extern int iter; //number of iterations

struct RNN //recurrent neural network
{
    int in; //number of inputs
    int hide; //number of feedback inputs
    matrix link; //connection matrix. Element [1][2] = 1, means net output 1 connected to feedback input 2
    int layers_amount; //number of layers (hide+output)
    matrix layers; //number of neurons in each layer
};

int RNN_par_amount(RNN net); //counting the number of RNN parameters

matrix RNN_par_init(RNN net); //initialization of RNN parameters

float norm_rand(float m, float q); //normal distribution

matrix giper_tan(const matrix& m); //activation function hyperbolic tangent

matrix ident(const matrix& m); //linear activation function

matrix RNN_calc(const matrix& IN, RNN net, const matrix& par, matrix (*func[])(const matrix&)); //get the output of RNN

matrix NN_calc(const matrix& IN, const matrix& layers, const matrix& par, matrix (*func[])(const matrix&)); //get the output of perceptron

matrix levenberg(matrix(*)(const matrix&, const matrix&, int),
                 const matrix& IN_exp, const matrix& OUT_exp,
                 const matrix& par0); //Levenberg-Marquardt algorithm

float mse(const matrix& a, const matrix& b); //mean square error betwen two data sets

matrix jacob(matrix(*)(const matrix&, const matrix&, int),
                 const matrix& IN_exp, const matrix& OUT_exp,
                 const matrix& par0); //Jacob matrix

#endif // NN_H
