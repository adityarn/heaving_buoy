class Linalg
{

    private:

// Function that returns the cross product
    // of r[], which is the positional vector from CG to control point
    // at x[j/n], y[j], z[j]
    double r_cross_n(int i, int j, double *cross)
    {
        double r[3];

        r[0] = (x[j/n] - cgx);

        r[1] = (y[j] - cgy);

        r[2] = (z[j] - cgz);

        if(i == 3)
            cross[0] = r[1]*nz[j] - r[2]*ny[j];

        if(i == 4)
            cross[1] = r[2]*nx[j/n] - r[0]*nz[j];

        if(i == 5)
            cross[2] = r[0]*ny[j] - r[1]*nx[j/n];

        return crossproduct;

    }

};
