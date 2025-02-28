#include "nn.h"



float h = 1e-4; //step to calculate derivatives
float u = 1; //regularization parameter
float u_max = 1e10; //max value of 'u'
float u_min = 1e-8; //min value of 'u'
float u_inc = 10; //coefficient to increase 'u'
float u_dec = 10; //coefficient to decrease 'u'
float err_min = 3e-4; //error when stop
int iter = 100; //number of iterations



int RNN_par_amount(RNN net)
{
    int out = 0;
    int a = net.in + net.hide; //number of inputs

    int i;
    for(i=0; i<net.layers_amount; i++)
    {
        out = out + (a + 1) * net.layers.read_data(0, i); //+weights to the next layer and biases for the next layer
        a = net.layers.read_data(0, i);
    }

    return out;
}



matrix RNN_par_init(RNN net)
{
    matrix par(1, RNN_par_amount(net));

    int a = net.in + net.hide; //number of inputs
    int b = 0;

    int i;
    int ii;
    int amount;
    for(i=0; i<net.layers_amount; i++) //for each layer
    {
        amount = a * net.layers.read_data(0, i); //number of weights

        float sigma; //standart deviation
        if(i<net.layers_amount-1)
        {
            sigma = 2 / (net.layers.read_data(0, i) + net.layers.read_data(0, i+1));
        }else
        {
            sigma = 1 / net.layers.read_data(0, i);
        }

        for(ii=0; ii<amount; ii++) //writing of weights with normal distribution
        {
            par.write_data(0, b+ii, norm_rand(0, sigma));
        }
        a = net.layers.read_data(0, i);
        b = b + amount + a;
    }

    return par;
}



float norm_rand(float m, float q) //normal distribution
{
    std::random_device generator; //true random number generator
    std::normal_distribution<float> distribution(m, q);

    return distribution(generator);
}


matrix giper_tan(const matrix& m)
{
    matrix A{m};

    int i;
    for(i=0; i<m.rows*m.columns; i++)
    {
        A.data[i] = tanh(m.data[i]);
    }

    return A;
}



matrix ident(const matrix& m)
{
    matrix A{m};

    return A;
}



matrix RNN_calc(const matrix& IN0, RNN net, const matrix& par0, matrix (*func[])(const matrix&))
{
    matrix IN = IN0;
    matrix par = par0;

    int t = IN.columns; //number of calculation cycles, time equivalent
    matrix OUT(net.layers.read_data(0, net.layers_amount-1), t); //output data
    matrix IN_for_perceptron(net.in+net.hide, 1); //input for perceptron

    //creating an array of layers for a perceptron, with all network input (data+feedback)
    int inputs[] = {net.in+net.hide};
    matrix A(1, 1, inputs);
    matrix layers = A.merg_by_columns(net.layers);

    int i;
    int ii;
    int iii;
    int a;
    for(i=0; i<t; i++) //calculation for each cycle
    {
        for(ii=0; ii<net.in; ii++) //write input data
        {
            IN_for_perceptron.write_data(ii, 0, IN.read_data(ii, i));
        }

        for(ii=0; ii<net.hide; ii++) //write feedback data
        {
            a = ii + net.in; //offset from beginning

            if(i>0) //if there is data for feedback, the feedback delay is 1 cycle
            {
                for(iii=0; iii<net.layers.read_data(0, net.layers_amount-1); iii++) //for each output
                {
                    if(net.link.read_data(iii, ii) == 1) //if there is feedback between the iii output and the ii feedback input
                    {
                        IN_for_perceptron.write_data(a, 0, OUT.read_data(iii, i-1)); //write the previous output state
                    }
                }
            }
        }

        OUT = OUT.insert_columns(i, NN_calc(IN_for_perceptron, layers, par, func));
    }

    return OUT;
}



matrix NN_calc(const matrix& IN0, const matrix& layers0, const matrix& par0, matrix (*func[])(const matrix&))
{
    matrix IN = IN0;
    matrix layers = layers0;
    matrix par = par0;

    int i;
    int ii;
    int iii;
    matrix w; //matrix of weights
    matrix b; //matrix of biases
    int adr = 0; //addresing by parameter array
    for(i=1; i<layers.columns; i++) //for each layer
    {
        w = matrix(layers.read_data(0, i), layers.read_data(0, i-1));
        for(ii=0; ii<layers.read_data(0, i); ii++) //write weights into the matrix
        {
            for(iii=0; iii<layers.read_data(0, i-1); iii++)
            {
                w.write_data(ii, iii, par.read_data(0, adr));
                adr++;
            }
        }

        b = matrix(layers.read_data(0, i), 1);
        for(ii=0; ii<layers.read_data(0, i); ii++) //write biases into the matrix
        {
            b.write_data(i, 0, par.read_data(0, adr));
            adr++;
        }

        matrix z = w * IN;
        IN = z + b;

        IN = func[i-1](IN);
    }

    return IN;
}



matrix levenberg(matrix(*model)(const matrix&, const matrix&, int),
                 const matrix& IN_exp, const matrix& OUT_exp,
                 const matrix& par0)
{
    matrix IN = IN_exp;
    matrix OUT = OUT_exp;
    matrix par = par0;
    matrix old_par;

    int i;
    float err0;
    float err;
    matrix r;
    matrix model_OUT;
    matrix J;
    matrix jj;
    matrix je;
    matrix I;
    matrix dp;
    for(i=0; i<iter; i++)
    {
        old_par = par;

        model_OUT = model(IN, par, 1); //get the current output of the model
        err0 = mse(OUT, model_OUT); //get the current MSE
        r = OUT - model_OUT; //deviations between data
        r = r.expand(); //transform into the row
        r = r.transpose(); //transform into the column

        J = jacob(model, IN, OUT, par); //Jacob matrix

        while(1)
        {
            jj = J.transpose() * J;
            je = J.transpose() * r;
            I = matrix(par.columns, par.columns).eye();
            dp = !(jj + I * u) * je;
            par = par - dp.transpose();

            model_OUT = model(IN, par, 1);
            err = mse(OUT, model_OUT);

            if(err<err0)
            {
                u = u / u_dec;
                if(u<u_min) u = u_min;
                break;
            }else
            {
                par = old_par;
                u = u * u_inc;
            }

            if(u>=u_max) return par;
        }

        err = mse(OUT, model(IN, par, 1));
        std::cout << "Iteration: " << i+1 << ", MSE: " << err << ", u: " << u << std::endl;

        if(err<err_min) break;
    }

    return par;
}



float mse(const matrix& a, const matrix& b)
{
    float err = 0;

    if(a.rows==b.rows && a.columns==b.columns)
    {
        int i;
        int ii;
        for(i=0; i<a.rows; i++)
        {
            for(ii=0; ii<a.columns; ii++)
            {
                err = err + powf(a.data[i*a.columns+ii] - b.data[i*a.columns+ii], 2);
            }
        }

        err = err / (float)(a.columns*a.rows-1);
        err = sqrtf(err);
    }else
    {
        std::cout << "mse(): Matrices of different size" << std::endl;
    }

    return err;
}



matrix jacob(matrix(*model)(const matrix&, const matrix&, int),
                 const matrix& IN_exp, const matrix& OUT_exp,
                 const matrix& par0)
{
    matrix OUT = OUT_exp;
    matrix IN = IN_exp;
    matrix par = par0;

    matrix J(OUT.expand().columns, par.columns);

    int i;
    matrix a;
    matrix r1;
    matrix r2;
    for(i=0; i<par.columns; i++)
    {
        par.data[i] += h;
        a = model(IN, par, 1);
        r1 = OUT - a;

        par.data[i] -= 2*h;
        a = model(IN, par, 1);
        r2 = OUT - a;

        a = r1 - r2;
        a = a.expand().transpose();
        a = a * 0.5 * (float)(1/h);

        par.data[i] += h;
        J = J.insert_columns(i, a);
    }

    return J;
}
