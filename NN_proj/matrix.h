#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>



class matrix
{
public:
    unsigned int rows;
    unsigned int columns;
    float* data;

    matrix(); //constructor without parameters
    matrix(unsigned int rows_in, unsigned int columns_in); //constructor with parameters
    matrix(unsigned int rows_in, unsigned int columns_in, float* data_in); //constructor with parameters
    matrix(unsigned int rows_in, unsigned int columns_in, int* data_in); //constructor with parameters
    matrix(const matrix& m); //copy constructor

    void print(); //ptint matrix
    void file_read(char* file_name); //read data (not raws and columns amont) from file
    void file_write(char* file_name); //write data in file
    void write_data(int i, int j, float d); //write data to [i][j], indexing from 0
    float read_data(int i, int j); //read data from [i][j], indexing from 0

    matrix& operator=(const matrix& m); //operator overloadng '=' matrix = matrix
    matrix operator+(const matrix& m); //operator overloadng '+' matrix + matrix
    matrix operator+(float n); //operator overloadng '+' matrix + number
    matrix operator*(const matrix& m); //operator overloadng '*' matrix * matrix
    matrix operator*(float n); //operator overloadng '*' matrix * number
    matrix operator-(const matrix& m); //operator overloadng '-' matrix - matrix
    matrix operator-(float n); //operator overloadng '-' matrix - number

    matrix transpose(); //matrix transpose
    matrix minor(int row, int column); //minor matrix of element in 'row' and 'column', indexing from 0
    double deter(); //determinant of matrix
    matrix operator!(); //inverse matrix
    matrix eye(); //identity matrix

    matrix merg_by_columns(const matrix& a); //merging of two matrices with preserving the number of rows
    matrix merg_by_rows(const matrix& a); //merging of two matrices with preserving the number of columns
    matrix shift(int* shifts); //data shift
    matrix cut_rows(int start, int stop); //data split into rows from start to finish, indexing from 0
    matrix cut_columns(int start, int stop); //data split into columns from start to finish, indexing from 0
    matrix insert_columns(int column, const matrix& a); //data rewriting with matrix 'a' from 'column', indexing from 0
    matrix expand(); //unfolds a matrix into a row
    float max(); //find max element
    float min(); //find min element

    ~matrix(); //destructor
};



#endif // MATRIX_H
