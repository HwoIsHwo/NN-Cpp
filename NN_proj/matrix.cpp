#include "matrix.h"



matrix::matrix()
{
    rows = 0;
    columns = 0;
    data = new float[0];
}



matrix::matrix(unsigned int rows_in, unsigned int columns_in) //constructor with parameters
{
    rows = rows_in;
    columns = columns_in;
    data = new float[rows*columns]{};
}



matrix::matrix(unsigned int rows_in, unsigned int columns_in, float* data_in)//constructor with parameters
{
    rows = rows_in;
    columns = columns_in;
    data = new float[rows*columns];

    int i;
    for(i=0; i<rows*columns; i++)
    {
        data[i] = data_in[i];
    }
}



matrix::matrix(unsigned int rows_in, unsigned int columns_in, int* data_in)//constructor with parameters
{
    rows = rows_in;
    columns = columns_in;
    data = new float[rows*columns];

    int i;
    for(i=0; i<rows*columns; i++)
    {
        data[i] = data_in[i];
    }
}



matrix::matrix(const matrix& m) //copy constructor
{
    rows = m.rows;
    columns = m.columns;
    data = new float[rows*columns];

    int i;
    for(i=0; i<rows*columns; i++)
    {
        data[i] = m.data[i];
    }
}



void matrix::print() //matrix output
{
    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<this->columns; ii++)
        {
            std::cout << data[i*columns+ii] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



void matrix::file_read(char* filename) //read data (not raws and columns amont) from file
{
    std::ifstream file;
    file.open(filename);
    if(file.is_open()!=1) //if file isn't opened
    {
        std::cout << "file_read(): File wasn't opened" << std::endl;
        return;
    }

    float number;
    int counter = 0;
    while(file>>number) counter++; //amount of data in the file
    file.close();

    *this = matrix(1, counter);

    file.open(filename);
    if(file.is_open()!=1) //if file isn't opened
    {
        std::cout << "file_read(): File wasn't opened" << std::endl;
        return;
    }

    int i;
    for(i=0; i<rows*columns; i++) //trying to read each element
    {
        if(file.eof()) break; // go out if it's the end of file
        file >> data[i];
    }

    file.close();
}



void matrix::file_write(char* filename) //write data in file
{
    std::ofstream file;
    file.open(filename);

    if(file.is_open()) //if file is opened
    {
        int i;
        for(i=0; i<rows*columns; i++) //write each element
        {
            file << data[i] << std::endl;
        }
    }else
    {
        std::cout << "file_write(): File wasn't opened" << std::endl;
    }

    file.close();
}



void matrix::write_data(int i, int j, float d)
{
    data[i*columns+j] = d;
}



float matrix::read_data(int i, int j)
{
    return data[i*columns+j];
}



matrix& matrix::operator=(const matrix& m) //operator overloadng '=' matrix = matrix
{
    if(&m != this) //check for writing in itself
    {
        delete [] data; //delete obj
        data = nullptr;

        rows = m.rows;
        columns = m.columns;
        data = new float[rows*columns];

        int i;
        for(i=0; i<rows*columns; i++) //copy data
        {
            data[i] = m.data[i];
        }
    }else
    {
        std::cout << "operator '=': Rewriting matrix by itself" << std::endl;
    }

    return *this;
}



matrix matrix::operator+(const matrix& m) //operator overloadng '+' matrix + matrix
{
    matrix A{*this}; //output matrix

    if(rows==m.rows && columns==m.columns) //check for dimensional conformity
    {
        int i;
        for(i=0; i<rows*columns; i++)
        {
            A.data[i] += m.data[i];
        }
    }else
    {
        std::cout << "operator 'm+m': Dimensions aren't equal" << std::endl;
    }

    return A;
}



matrix matrix::operator+(float n) //operator overloadng '+' matrix + number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++)
    {
        A.data[i] += n;
    }

    return A;
}



matrix matrix::operator*(const matrix& m) //operator overloadng '*' matrix * matrix
{
    matrix A(rows, m.columns); //output matrix

    if(columns == m.rows) //check for dimensional conformity
    {
        int i;
        int ii;
        int iii;
        float element;
        for(i=0; i<A.rows; i++) //for each row of A matrix
        {
            for(ii=0; ii<A.columns; ii++) //for each element in row
            {
                element = 0;
                for(iii=0; iii<columns; iii++) //calculation of product i-row of left matrix and ii-column of right matrix
                {
                    element = element + data[i*columns+iii] * m.data[iii*m.columns+ii];
                }
                A.data[i*A.columns+ii] = element;
            }
        }
    }else
    {
        std::cout << "operator 'm*m': Columns amount of left matrix != Rows amount of right matrix" << std::endl;
    }

    return A;
}



matrix matrix::operator*(float n) //operator overloadng '*' matrix * number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++) //copy data
    {
        A.data[i] *= n;
    }

    return A;
}



matrix matrix::operator-(const matrix& m) //operator overloadng '-' matrix - matrix
{
    matrix A{*this}; //output matrix

    if(rows==m.rows && columns==m.columns) //check for dimensional conformity
    {
        int i;
        for(i=0; i<rows*columns; i++)
        {
            A.data[i] -= m.data[i];
        }
    }else
    {
        std::cout << "operator 'm-m': Dimensions aren't equal" << std::endl;
    }

    return A;
}



matrix matrix::operator-(float n) //operator overloadng '-' matrix - number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++)
    {
        A.data[i] -= n;
    }

    return A;
}



matrix matrix::transpose() //matrix transpose
{
    matrix A(columns, rows); //output matrix

    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<columns; ii++)
        {
            A.data[ii*rows+i] = data[i*columns+ii];
        }
    }

    return A;
}



matrix matrix::minor(int row, int column) //minor matrix of element in 'row' and 'column'
{
    if(rows==columns) //matrix is square
    {
        matrix A(rows-1, columns-1);

        int i;
        int ii;
        int cnt = 0;
        for(i=0; i<rows; i++)
        {
            if(i==row) continue;

            for(ii=0; ii<columns; ii++)
            {
                if(ii==column) continue;

                A.data[cnt] = data[i*columns+ii];
                cnt++;
            }
        }

        return A;
    }else
    {
        std::cout << "minor(): Matrix isn't square" << std::endl;
        matrix A{*this};

        return A;
    }
}



double matrix::deter() //determinant of matrix
{
    double D = 1;
    matrix A{*this};

    if(rows==columns) //matrix is square
    {
        int i;
        int ii;
        int iii;
        float k;
        for(i=0; i<rows; i++)
        {
            for(ii=1+i; ii<rows; ii++)
            {
                k = A.data[ii*columns+i] / A.data[i*columns+i];

                for(iii=0+i; iii<columns; iii++)
                {
                    A.data[ii*columns+iii] -= A.data[i*columns+iii] * k;
                }
            }
        }

        for(i=0; i<rows; i++)
        {
            D = D * A.data[i*columns+i];
        }
    }else
    {
        std::cout << "deter(): Matrix isn't square" << std::endl;
    }

    return D;
}



matrix matrix::operator!() //inverse matrix
{
    matrix A{*this}; //output matrix

    if(rows==columns) //matrix is square
    {
        //double d = deter();
        double d = 1;
        if(d!=0)
        {
            A = A.merg_by_columns(A.eye()); //add the identity matrix

            int i;
            int ii;
            int iii;
            float k;
            for(i=0; i<rows; i++)
            {
                k = A.data[i*A.columns+i];
                for(ii=i; ii<A.columns; ii++)
                {
                    A.data[i*A.columns+ii] /= k;
                }

                for(ii=0; ii<rows; ii++)
                {
                    if(ii==i) continue;
                    k = A.data[ii*A.columns+i];
                    for(iii=0; iii<A.columns; iii++)
                    {
                        A.data[ii*A.columns+iii] -= A.data[i*A.columns+iii] * k;
                    }
                }
            }

            A = A.cut_columns(columns, 2*columns-1);
        }else
        {
            std::cout << "operator '!': Determinant = 0, it's impossible to get iverse matrix" << std::endl;
        }
    }else
    {
        std::cout << "operator '!': Matrix isn't square" << std::endl;
    }

    return A;
}



matrix matrix::eye()
{
    matrix A{*this};

    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<columns; ii++)
        {
            if(i==ii)
            {
                A.data[i*columns+ii] = 1;
            }else
            {
                A.data[i*columns+ii] = 0;
            }
        }
    }

    return A;
}



matrix matrix::merg_by_columns(const matrix& a)
{
    matrix A(rows, a.columns+columns);
    if(a.rows==rows) //the number of matrix rows is equal
    {
        int i;
        int ii;
        for(i=0; i<rows; i++) //copy data from first matrix
        {
            for(ii=0; ii<columns; ii++)
            {
                A.data[i*A.columns+ii] = data[i*columns+ii];
            }
        }

        for(i=0; i<rows; i++) //copy data from second matrix
        {
            for(ii=0; ii<a.columns; ii++)
            {
                A.data[i*A.columns+ii+columns] = a.data[i*a.columns+ii];
            }
        }
    }else
    {
        std::cout << "merg_by_columns(): The number of matrix rows isn't equal" << std::endl;
    }

    return A;
}



matrix matrix::merg_by_rows(const matrix& a)
{
    matrix A(a.rows+rows, columns);
    if(a.columns==columns) //the number of matrix columns is equal
    {
        int i;
        int ii;
        for(i=0; i<rows; i++) //copy data from first matrix
        {
            for(ii=0; ii<columns; ii++)
            {
                A.data[i*A.columns+ii] = data[i*columns+ii];
            }
        }

        for(i=0; i<a.rows; i++) //copy data from second matrix
        {
            for(ii=0; ii<a.columns; ii++)
            {
                A.data[(i+rows)*A.columns+ii] = a.data[i*a.columns+ii];
            }
        }
    }else
    {
        std::cout << "merg_by_rows(): The number of matrix columns isn't equal" << std::endl;
    }

    return A;
}



matrix matrix::shift(int* shifts)
{
    matrix A(rows, columns);

    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=shifts[i]; ii<columns; ii++)
        {
            A.data[i*columns+ii] = data[i*columns+ii-shifts[i]];
        }
    }

    return A;
}



matrix matrix::cut_rows(int start, int stop)
{
    matrix A(stop-start+1, columns);

    if(start<=stop && start>=0 && stop<rows)
    {
        int i;
        int ii;
        for(i=start; i<=stop; i++) //for rows from start to stop
        {
            for(ii=0; ii<columns; ii++)
            {
                A.data[(i-start)*columns+ii] = data[i*columns+ii];
            }
        }
    }else
    {
        std::cout << "cut_rows(): Incorrect borders" << std::endl;
    }

    return A;
}



matrix matrix::cut_columns(int start, int stop)
{
    matrix A(rows, stop-start+1);

    if(start<=stop && start>=0 && stop<columns)
    {
        int i;
        int ii;
        for(i=0; i<rows; i++) //for rows from start to stop
        {
            for(ii=start; ii<=stop; ii++)
            {
                A.data[i*A.columns+ii-start] = data[i*columns+ii];
            }
        }
    }else
    {
        std::cout << "cut_columns(): Incorrect borders" << std::endl;
    }

    return A;
}



matrix matrix::insert_columns(int column, const matrix &a)
{
    int i;
    int ii;
    matrix A{*this};

    if(rows==a.rows && (columns-column)>=a.columns)
    {
        for(i=0; i<a.rows; i++)
        {
            for(ii=0; ii<a.columns; ii++)
            {
                A.data[i*columns+ii+column] = a.data[i*a.columns+ii];
            }
        }
    }else
    {
        std::cout << "insert_columns(): Number of rows isn't equal or array out of bounds" << std::endl;
    }

    return A;
}



matrix matrix::expand()
{
    matrix A{*this};

    A.columns = A.columns * A.rows;
    A.rows = 1;

    return A;
}



float matrix::max()
{
    float out = 1e-8;

    int i;
    for(i=0; i<rows*columns; i++)
    {
        if(data[i]>out) out = data[i];
    }

    return out;
}



float matrix::min()
{
    float out = 1e8;

    int i;
    for(i=0; i<rows*columns; i++)
    {
        if(data[i]<out) out = data[i];
    }

    return out;
}



matrix::~matrix() //destructor
{
    delete [] data;
    data = nullptr;
}
