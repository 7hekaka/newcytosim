// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University

#ifndef ACCUMULATOR_H
#define ACCUMULATOR_H

/// Helper class to calculate moments of a cloud of points
/**
This calculates the mean and the variance of a set of 3D points
 @todo update to non-naive method that also work with large samples
 https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
 */
class Accumulator
{
public:
    double sum;
    double avg[3];
    double var[9];
    
    /// set all accumulators to zero
    void reset()
    {
        sum = 0;
        for ( int i = 0; i < 3; ++i ) avg[i] = 0;
        for ( int i = 0; i < 9; ++i ) var[i] = 0;
    }
    
    Accumulator()
    {
        reset();
    }

    /// add `p` with weight 'w'
    void add(real w, Vector const& p)
    {
        sum += w;
        avg[0] += w * p.XX;
        var[0] += w * p.XX * p.XX;
#if ( DIM > 1 )
        avg[1] += w * p.YY;
        var[1] += w * p.YY * p.XX;
        var[4] += w * p.YY * p.YY;
#endif
#if ( DIM > 2 )
        avg[2] += w * p.ZZ;
        var[2] += w * p.ZZ * p.XX;
        var[5] += w * p.ZZ * p.YY;
        var[8] += w * p.ZZ * p.ZZ;
#endif
    }
    
    /// add `p` with weight 1.0
    void add(Vector const& p)
    {
        sum += 1;
        avg[0] += p.XX;
        var[0] += p.XX * p.XX;
#if ( DIM > 1 )
        avg[1] += p.YY;
        var[1] += p.YY * p.XX;
        var[4] += p.YY * p.YY;
#endif
#if ( DIM > 2 )
        avg[2] += p.ZZ;
        var[2] += p.ZZ * p.XX;
        var[5] += p.ZZ * p.YY;
        var[8] += p.ZZ * p.ZZ;
#endif
    }
    
    /// transform the second-order acumulators into variance
    void subtract_mean()
    {
        const double scale = 1.0 / ( sum - 1 );
        avg[0] /= sum;
        var[0] = var[0]*scale - avg[0] * avg[0];
#if ( DIM > 1 )
        avg[1] /= sum;
        var[1] = var[1]*scale - avg[1] * avg[0];
        var[4] = var[4]*scale - avg[1] * avg[1];
#endif
#if ( DIM > 2 )
        avg[2] /= sum;
        var[2] = var[2]*scale - avg[2] * avg[0];
        var[5] = var[5]*scale - avg[2] * avg[1];
        var[8] = var[8]*scale - avg[2] * avg[2];
#endif
    }
    
    real total_length() const
    {
        return sum;
    }
    
    real total_variance() const
    {
        return var[0] + var[4] + var[8];
    }
    
    void print_doc(std::ostream& os) const
    {
        os << COM << "cnt";
        os << SEP << "avgX" << SEP << "avgY" << SEP << "avgZ";
        os << SEP << "varX" << SEP << "varY" << SEP << "varZ";
        os << SEP << "var_sum";
    }
    
    void print(std::ostream& os, bool mode)
    {
        if ( mode )
            os << LIN << (int) sum;
        else
            os << SEP << sum;
        os << SEP << avg[0];
        os << SEP << avg[1];
        os << SEP << avg[2];
        os << SEP << var[0];
        os << SEP << var[4];
        os << SEP << var[8];
        os << SEP << var[0] + var[4] + var[8];
    }
};

#endif
