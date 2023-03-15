#!/usr/bin/env python3
# `go_sim_lib.py` is a miniature library to run Cytosim.
#  It is not executable directly, and instead it is used by go_sim.py
#  to create a directory, copy files, move directories, etc.
#
# Copyright F. Nedelec 2007--2023 with S. Dmitrieff 2019


try:
    import os, sys, shutil, subprocess
except ImportError:
    host = os.getenv('HOSTNAME', 'unknown')
    sys.stderr.write("go_sim_lib.py could not load python modules on %s\n" % host)
    sys.exit()

try:
    import exceptions
except:
    try:
        import builtins as exceptions
    except:
        host = os.getenv('HOSTNAME', 'unknown')
        sys.stderr.write("go_sim_lib.py could not load `exceptions` on %s\n" % host)
        sys.exit()

class Error( exceptions.Exception ):
    """go_sim.py exception class"""
    def __init__(self, value=None):
        self.value = value
    def __str__(self):
        return repr(self.value)

# name of the log file:
logfile_name = 'log.txt'

# default output for error messages:
err = sys.stderr
out = sys.stdout

#==========================  DIR/FILES HANDLING ==============================

def make_directory(root, n=0):
    """
    Create a new directory `root####`, using a 4-digit number >= n
    """
    if root[-1].isdigit():
        try:
            os.mkdir(root)
            return root
        except OSError:
            root = root + '-'
    while n < 10000:
        res = root + '%04i' % n
        try:
            os.mkdir(res)
            return res
        except OSError:
            pass #err.write("failed " + res)
        n += 1
    raise Error("failed to create new run directory on "+os.getenv('HOSTNAME', 'unknown'))


def make_run_directory(root, conf):
    """create a temporary directory starting by `root`"""
    import tempfile
    if 'SLURM_JOB_ID' in os.environ:
        try:
            # RDS directory on Cambridge's Research Computing Services
            path = os.path.dirname(conf)
            if path.endswith('todo'):
                path = path[:-4]
            return tempfile.mkdtemp('', root+'-', path)
        except:
            pass
        try:
            return tempfile.mkdtemp('', root, '/local')
        except:
            pass
        tmp = os.getenv('TMPDIR', '')
        if tmp:
            return tmp
    return make_directory(root)
    #return tempfile.mkdtemp('', root, '.')


def copy_recursive(src, dst):
    """recursively copy everything from src to dst"""
    if os.path.isfile(src):
        shutil.copy2(src, dst)
    elif os.path.isdir(src):
        try:
            os.mkdir(dst)
        except OSError:
            pass
        files = os.listdir(src)
        for f in files:
            s = os.path.join(src, f)
            d = os.path.join(dst, f)
            copy_recursive(s, d)


def park_directory(path, park, name):
    """Copy directory 'path' to park, under a similar name"""
    src = os.path.abspath(path)
    dst = os.path.join(park,name)
    if src == os.path.abspath(park):
        return src
    try:
        shutil.copytree(src, dst)
    except:
        dst = make_directory(dst)
        copy_recursive(src, dst)
    out.write("moving ( %s -> %s )" % (src, dst))
    from filecmp import dircmp
    dcmp = dircmp(src, dst)
    if dcmp.left_only or dcmp.diff_files:
        out.write(" ---> failed!" % path)
        err.write("go_sim_lib.py failed to copy '%s' verbatim\n" % path)
        return src
    else:
        shutil.rmtree(src)
    out.write(" : done\n")
    return dst


def copy_config(name, repeat):
    """
        make 'repeat' copies of the name.
    """
    res = []
    for x in range(repeat):
        res.extend([name])
    return res


def make_config(conf, repeat, script, dest):
    """
    Generate config files by running a python script,
    or simply repeat the name if ( repeat > 1 ) and script==''.
    """
    if script:
        code = script.rstrip('.py')
        module = {}
        try:
            module = __import__(code)
        except:
            import imp
            module = imp.load_source(code, script)
        if not module:
            raise Error("could not load python module `"+code+"'")
        # use module to generate a new config file:
        return module.parse(conf, {}, repeat, dest)
    else:
        return copy_config(conf, repeat)


#=======================  RUNNING THE SIMULATION  ==============================

def run_sim(exe, args):
    """
    Start executable in current directory, and wait for completion.
    Standard output is sent to `out.txt' and standard error to `err.txt'.
    The executable shall find its default configuration file.
    """
    outname = 'out.txt'
    errname = 'err.txt'
    outfile = open(outname, 'w')
    errfile = open(errname, 'w')
    # run simulation
    if not args:
        val = subprocess.call(exe, stdout=outfile, stderr=errfile)
    else:
        val = subprocess.call([exe]+args, stdout=outfile, stderr=errfile)
    outfile.close()
    errfile.close()
    # remove output files if empty:
    if os.path.isfile(outname) and not os.path.getsize(outname):
        os.remove(outname)
    if os.path.isfile(errname) and not os.path.getsize(errname):
        os.remove(errname)
    return val


def info_start(filename, exe, conf, args, pid):
    import time
    with open(filename, "w") as f:
        f.write("host      %s\n" % os.getenv('HOSTNAME', 'unknown'))
        f.write("user      %s\n" % os.getenv('USER', 'unknown'))
        f.write("wdir      %s\n" % os.getcwd())
        f.write("exec      %s\n" % exe)
        f.write("args      %s\n" % args)
        f.write("conf      %s\n" % conf)
        f.write("pid       %s\n" % pid)
        f.write("start     %s\n" % time.asctime())


def info_end(filename, val):
    import time
    with open(filename, "a") as f:
        f.write("status    %s\n" % val)
        f.write("stop      %s\n" % time.asctime())


def run(exe, conf, args, job_name, exe_input_name):
    """
    Run one simulation in a new sub directory and wait for completion.
    The config file 'conf' is copied to the subdirectory.
    Returns sub-directory in which `exe` was called.
    """
    cdir = os.getcwd()
    if not os.path.isfile(conf):
        raise Error("missing/unreadable config file")
    conf = os.path.abspath(conf);
    wdir = make_run_directory(job_name, conf)
    os.chmod(wdir, 504)
    os.chdir(wdir)
    shutil.copyfile(conf, exe_input_name)
    info_start(logfile_name, exe, conf, args, os.getpid())
    val = run_sim(exe, args)
    info_end(logfile_name, val)
    os.chdir(cdir)
    return (val, wdir)


def start(exe, conf, args, root, exe_input_name):
    """
    Start simulation in a new sub directory, and return immediately.
    The config file `conf` is copied to the sub-directory.
    """
    cdir = os.getcwd()
    if not os.path.isfile(conf):
        raise Error("missing/unreadable config file")
    conf = os.path.abspath(conf)
    wdir = make_directory(root)
    os.chdir(wdir)
    shutil.copyfile(conf, exe_input_name)
    if exe != 'none':
        outfile = open('out.txt', 'w')
        errfile = open('err.txt', 'w')
        #start simulation, but do not wait for completion:
        pid = subprocess.Popen(['nohup', exe]+args, stdout=outfile, stderr=errfile).pid
        info_start(logfile_name, exe, conf, args, pid)
    else:
        pid = 0
    os.chdir(cdir)
    return (pid, wdir)


