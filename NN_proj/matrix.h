#ifndef MATRIX_H
#define MATRIX_H


class matrix
{
public:
    unsigned int rows;
    unsigned int columns;
    float* data;

    matrix(unsigned int rows_in, unsigned int columns_in); //constructor with parameters
    matrix(unsigned int rows_in, unsigned int columns_in, float* data_in); //constructor with parameters
    matrix(const matrix &m); //copy constructor

    void print(); //ptint matrix
    void file_read(char* file_name); //read data (not raws and columns amont) from file
    void file_write(char* file_name); //write data in file
    void write_data(int i, int j, float d);
    float read_data(int i, int j);

    matrix& operator=(const matrix& m); //operator overloadng '=' matrix = matrix
    matrix operator+(const matrix& m); //operator overloadng '+' matrix + matrix
    matrix operator+(float n); //operator overloadng '+' matrix + number
    matrix operator*(const matrix& m); //operator overloadng '*' matrix * matrix
    matrix operator*(float n); //operator overloadng '*' matrix * number
    matrix operator-(const matrix& m); //operator overloadng '-' matrix - matrix
    matrix operator-(float n); //operator overloadng '-' matrix - number

    matrix transpose(); //matrix transpose
    matrix minor(int row, int column); //minor matrix of element in 'row' and 'column'
    float deter(); //determinant of matrix
    matrix operator!(); //inverse matrix
    void eye(); //identity matrix

    ~matrix(); //destructor
};

#endif // MATRIX_H
