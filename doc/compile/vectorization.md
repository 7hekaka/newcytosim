
# SIMD Vectorization

Many routines in cytosim are implemented using vectorized primitives to achieve significant speedup on processors that support this technology. It uses the SIMD instruction set, to perform operations on multiple scalar values in parallel.
For example 'AVX' allow to operate on vectors of 4 double precision floats.

There has been different generations of SIMD, of increased vector size:

    128 bits:  SSE3           ~2004
    256 bits:  AVX            ~2009
    256 bits:  AVX2 and FMA   ~2013
    512 bits:  AVX-512        ~2018

Typically, one can expect a 40% speed up of AVX compared to SSE3, and 
another speedup with the 'Fused Add-multiply' instructions of `FMA`.
The functions included in `AVX2` other than FMA are not used in cytosim.
To take advantage of this, you need to adjust the `makefile.inc`.

To check SIMD, you can run `test_matrix` to compare timing and checksums:

    Matrix test and timing code --- real 8 --- 4.2.1 Compatible Apple LLVM 11.0.0 (clang-1100.0.33.17)
    ------ 2D size 4080  filled 1.6 % :
    SMS1e 992208             set     18.37  muladd     28.53  alt     28.33  check     -2077.026722     -2077.026722 
    SMS2e 992208             set     15.74  muladd     29.19  alt     29.45  check     -2077.026722     -2077.026722 
    SMSBe 4*248562           set      7.67  muladd     16.94  alt     18.94  check     -2077.026722     -2077.026722 
    SMSBD 4*250602           set      8.01  muladd     17.81  alt     18.85  check     -2077.026722     -2077.026722 
    SMB 4*495084             set      7.35  muladd     23.96  alt     24.38  check     -2077.026722     -2077.026722 

SMS1 is "Sparse Matrix Symmetric1". In Cytosim, Meca is normally using SMS1 and SMSB only.
Here all the matrices were compiled with SSE3, which is indicated with the 'e' added to the name.
It would be an 'x' when AVX is used:

    Matrix test and timing code --- real 8 --- 4.2.1 Compatible Apple LLVM 11.0.0 (clang-1100.0.33.17)
    ------ 2D size 4080  filled 1.6 % :
    SMS1x 991148             set     15.20  muladd     31.33  alt     30.96  check     -1125.230053     -1125.230053 
    SMS2e 991148             set     14.40  muladd     30.03  alt     30.17  check     -1125.230053     -1125.230053 
    SMSBx +4*248297          set      8.21  muladd     15.01  alt     19.28  check     -1125.230053     -1125.230053 
    SMSBDx +4*250337         set      7.66  muladd     14.08  alt     18.93  check     -1125.230053     -1125.230053 
    SMBx +4*494554           set      7.13  muladd     22.49  alt     24.45  check     -1125.230053     -1125.230053 

The compiler might automatically promote SSE3 instructions into AVX ones.

# Compilation

Edit ‘makefile.inc’ to set the correct option for your compiler.
For example, to use `AVX` and `FMA` with `gcc`:

    -mavx -mfma

And for SSE3:

    -msse3

With Intel compiler (icpc):

    # Intel advanced instruction sets:
    # '-xHost' to optimize for host machine
    # '-xAVX' for AVX
    # '-march=core-avx2' for Intel core i7 (ca. 2015)

# SLURM submission

In `submit_slurm.py`, make sure this option is set in `sub()`

    def sub(exe):
        """return command that will submit one job"""
        # specify memory, shell, minimum number of cores and queue
        cmd  = [subcmd, '--nodes=1', '--ntasks=1']
        ...
  	   	 
and also for any job arrays in `array()`:

    def array(jobcnt):
        """return command that will submit a job-array"""
        # define parameters directly in the script:
        cmd  = ['#SBATCH --nodes=1']
        cmd += ['#SBATCH --ntasks=1']
        ...

# Troubleshooting

The job should run. Otherwise, you will get 'Illegal instruction' if the instructions are not supported by the CPU on which the program is running. In this case, the program will terminate with UNIX signal `SIGILL` of value 4.


FJN 27/04/2018
