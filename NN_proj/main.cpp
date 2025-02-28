#include <iostream>
#include <termios.h>
#include <sys/ioctl.h>
#include "matrix.h"
#include "nn.h"



matrix model(const matrix& IN, const matrix& par, int flag) //model function, flag=0 - init of parameters, flag=1 - get function output
{
    matrix out;

    RNN net1; //see nn.h
    net1.in = 1;
    net1.hide = 1;
    net1.link = matrix(1, 1).eye();
    net1.layers_amount = 2;
    float l1[]{10, 1};
    net1.layers = matrix(1, 2, l1);

    if(flag==0)
    {
        out = RNN_par_init(net1);
    }else
    {
        matrix A = IN; //data shift
        int shifts[] = {1};
        A.shift(shifts);

        matrix (*func[])(const matrix&) = {giper_tan, ident}; //activation function for each layer

        out = RNN_calc(A, net1, par, func);
    }

    return out;
}



void plot(const matrix& time, const matrix& y1, const matrix& y2)
{
    struct winsize ws;
    int rc;
    rc = ioctl(0, TIOCGWINSZ, &ws); //get the console size

    int row = 24 - 1;
    int col = 80;
    if(rc>=0)
    {
        row = ws.ws_row - 1; //-1 empty row
        col = ws.ws_col;
    }
    char* image = new char[row*col]{};
    matrix x = time;
    matrix y = y1;
    matrix y0 = y2;

    //axes
    int i;
    int ii;
    for(i=0; i<row*col; i++)
    {
        image[i] = ' ';
    }
    for(i=0; i<row; i++)
    {
        image[i*col] = '|';
    }
    image[(row-1)*col] = 'L';
    for(i=1; i<col; i++)
    {
        image[(row-1)*col+i] = '_';
    }

    //averaging data 1
    int window = x.columns / (col-1); //data averaging window
    matrix compressed_data_y(1, col-1); //matrix for compressed data (after averaging)
    float aver = 0;
    for(i=0; i<col-1; i++)
    {
        for(ii=0; ii<window; ii++)
        {
            aver += y.read_data(0, i*window+ii);
        }
        aver /= (float)window;
        compressed_data_y.write_data(0, i, aver);
        aver = 0;
    }

    //write data 1 to image
    float h = compressed_data_y.max()*1.1/(row-1); //step of the vertical scale
    for(i=0; i<col-1; i++)
    {
        ii = (int)(compressed_data_y.read_data(0, i) / h);
        image[(row-1-ii)*col+i] = '*';
    }

    //averaging data 2
    aver = 0;
    for(i=0; i<col-1; i++)
    {
        for(ii=0; ii<window; ii++)
        {
            aver += y0.read_data(0, i*window+ii);
        }
        aver /= (float)window;
        compressed_data_y.write_data(0, i, aver);
        aver = 0;
    }

    //write data 2 to image
    for(i=0; i<col-1; i++)
    {
        ii = (int)(compressed_data_y.read_data(0, i) / h);
        image[(row-1-ii)*col+i] = '0';
    }

    //print
    for(i=0; i<row; i++)
    {
        for(ii=0; ii<col; ii++)
        {
            std::cout << image[i*col+ii];
        }
        std::cout << std::endl;
    }

    //getchar();

    delete [] image;
}



int main()
{
    char input_data[] = {"Data/X.txt"};
    char parameters[] = {"Data/par.txt"};
    char output_data[] = {"Data/Y.txt"};
    char time[] = {"Data/t.txt"};
    int n = 50;
    iter = 100;

    //read input data
    matrix IN;
    IN.file_read(input_data);
    IN = IN.cut_columns(0, (int)IN.columns/n);

    //read output data
    matrix OUT;
    OUT.file_read(output_data);
    OUT = OUT.cut_columns(0, (int)OUT.columns/n);

    //read time
    matrix t;
    t.file_read(time);
    t = t.cut_columns(0, (int)t.columns/n);

    //read params
    matrix par;
    par.file_read(parameters);

    /*//init params
    matrix par = model({}, {}, 0);
    par.file_write(parameters);*/

    //learning
    par = levenberg(model, IN, OUT, par);
    par.file_write(parameters);
    matrix model_OUT = model(IN, par, 1);
    plot(t, model_OUT, OUT);

    /*int i = 0;
    while(1)
    {
        if(i==100000000)
        {
            plot(t, model_OUT);
            i = 0;
        }else
        {
            i++;
        }
    }*/

    return 0;
}
