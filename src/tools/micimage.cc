// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 Generates pseudo light-microscope images, from coordinate files
 F. Nedelec, 6 March 2013
 */

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdint.h>

#include "random.h"
#include "glossary.h"
#include "filepath.h"
#include "save_image_gl.h"


size_t   hits = 1000;
size_t dim[2] = { 1024, 1024 };
real      pix = 6.5;
real      mag = 100;
real      res = pix / mag;
real   foc[3] = { 0, 0, 0 };
real      psf = 0.22;
real      bck = 0;
real      brt = 1;


//------------------------------------------------------------------------------

/** Set array with background noise */
void clear(uint16_t array[], size_t array_size)
{
    if ( bck > 0 )
    {
        if ( bck < 100 )
        {
            for ( size_t i = 0; i < array_size; ++i )
                array[i] = RNG.poisson(bck);
        }
        else
        {
            // we approximate the sum of Poisson with a Gaussian
            real * g = new_real(array_size);
            RNG.gauss_set(g, array_size);
            
            for ( size_t u = 0; u < array_size; ++u )
            {
                if ( g[u] > -1.0 )
                    array[u] = bck * ( 1.0 + g[u] );
                else
                    array[u] = 0;
            }
            
            free_real(g);
        }
    }
    else
    {
        for ( size_t i = 0; i < array_size; ++i )
            array[i] = 0;
    }
}


void emitter(uint16_t array[], size_t array_size, real x, real y, real z)
{
    for ( size_t p = 0; p < hits; ++p )
    {
        real px, py;
        RNG.gauss_set(px, py, psf);
        px = ( x + px ) / res;
        py = ( y + py ) / res;
        if ( 0 < px  &&  0 < py )
        {
            int ix = (int) px;
            int iy = (int) py;
            if ( ix < dim[0] && iy < dim[1] )
                array[iy+dim[0]*ix] += 1;
        }
    }
}


size_t image(std::istream& is, FILE * file)
{
    // lower left corner of image:
    real cx = foc[0] - 0.5 * res * dim[0];
    real cy = foc[1] - 0.5 * res * dim[1];
    
    size_t array_size = dim[0] * dim[1];
    uint16_t * array = new uint16_t[array_size];
    
    clear(array, array_size);
    
    size_t cnt = 0;
    std::string line;
    while ( is.good() )
    {
        std::getline(is, line);
        
        if ( line.size() < 3 || line[0]=='%' )
            continue;

        std::istringstream ls(line);
        
        real x = 0, y = 0, z = 0;
        ls >> x;
        
        if ( ls.fail() )
            continue;

        ls >> y >> z;
        
        if ( RNG.test(brt) )
        {
            emitter(array, array_size, x-cx, y-cy, z);
            ++cnt;
        }
    }
    
    SaveImage::saveGrayPNG(file, array, dim[0], dim[1]);
    delete[] array;
    return cnt;
}


//------------------------------------------------------------------------------

void help(std::ostream& os)
{
    os << "Synopsis: \n";
    os << "   Generates pseudo light-microscope images, from coordinate files\n";
    os << "\n";
    os << "Syntax:\n";
    os << "\n";
    os << "   micimage INPUT-FILE [MORE-FILES] [options]\n";
    os << "\n";
    os << "   Input files should should be text files specifying the positions of emitters,\n";
    os << "   one per line, as ' X Y ' or ' X Y Z '\n";
    os << "   Coordinates can be specified as floating points values.\n";
    os << "   Empty lines and lines starting with '%' are ignored\n";
    os << "   Multiple files can be specified, and micimage will generate one\n";
    os << "   16bit gray-level PNG image file for each input file\n";
    os << "\n";
    os << "Method:\n";
    os << "\n";
    os << "    Each emitter listed in the file has the probability 'bright' to be active\n";
    os << "    An active emitter generates 'hits' increments on the detector array.\n";
    os << "    The hits are spread using a Gaussian with a standard deviation 'psf' around\n";
    os << "    the position of the emitter (a Gaussian is used to approximate the hairy-discs)\n";
    os << "    Pixels record the number of hits, to which a background noise is added\n";
    os << "\n";
    os << "    OPTION           DEFAULT     DESCRIPTION\n";
    os << "    dim=INT,INT      1024,1024   nb of pixels ( width, height )\n";
    os << "    focus=REAL,REAL  0,0         coordinates of the point seen at the center\n";
    os << "    back=REAL        0           level of Poisson noise in each pixel\n";
    os << "    pix=REAL         6.5         real size of a pixel on the detector\n";
    os << "    mag=REAL         100         objective magnification\n";
    os << "\n";
    os << "    bright=REAL      1           probability for an emitter to be active\n";
    os << "    psf=REAL         0.22        standard-deviation of hit spread ~ 0.5*wavelength/N.A.\n";
    os << "    hits=INT         1000        nb of detector hits produced by each active emitter\n";
    os << "\n";
    os << "Example:\n";
    os << "\n";
    os << "reportF fiber:speckles interval=1\n";
    os << "micimage report*.txt\n";
    os << "\n";
    os << "F. Nedelec, 6.03.2013\n";
    
    std::endl(os);
}

//------------------------------------------------------------------------------

std::string filename(std::string const& name)
{
    std::string res;
    std::string::size_type p = name.find('.');
    if ( p != std::string::npos )
        res = name.substr(0, p+1);
    else
        res = name;
    return res.append("png");
}


typedef std::vector<std::string> FileList;

int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    FileList files;
    Glossary arg;

    int n = 1;
    while ( n < argc )
    {
        const char * arg = argv[n];
        if ( memchr(arg, '=', strlen(arg)) )
            break;
        files.push_back(argv[n]);
        ++n;
    }
     
    if ( arg.read_strings(argc-n, argv+n) )
        return EXIT_FAILURE;

    arg.set(mag,    "mag");
    arg.set(pix,    "pix");
    arg.set(bck,    "back");
    arg.set(brt,    "bright");
    arg.set(hits,   "hits");
    arg.set(foc, 2, "focus");
    arg.set(dim, 2, "dim");
    arg.set(psf,    "psf");

    arg.print_warnings(std::cerr, 1, "\n");

    res = pix / mag;
    std::cerr << "A pixel covers a distance of " << res << " units" << '\n';

    RNG.seed();
    
    for ( std::string const& f : files )
    {
        std::string in = FilePath::file_part(f);
        std::ifstream is(in.c_str(), std::ifstream::in);
        
        if ( is.good() )
        {
            std::string on = filename(in);
            FILE * file = fopen(on.c_str(), "wb");

            if ( file )
            {
                if ( !ferror(file) )
                {
                    size_t cnt = image(is, file);
                    std::cerr << "created " << on << " from " << cnt << " active emitters" << '\n';
                }
                fclose(file);
            }
        }
        else
        {
            std::cerr << "Could not open file " << in << '\n';
        }
        
        is.close();
    }
}
