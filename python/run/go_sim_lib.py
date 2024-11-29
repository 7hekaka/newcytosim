#!/usr/bin/env python3
# `go_sim_lib.py` is a miniature library to run Cytosim.
#  It is not executable directly, and instead it is used by go_sim.py
#  to create a directory, copy files, move directories, etc.
#
# Copyright F. Nedelec 2007--2024 with S. Dmitrieff 2019


try:
    import os, sys, shutil
    from subprocess import Popen
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
        # RDS directory on Cambridge's Research Computing Services
        path = os.path.dirname(conf)
        if path.endswith('/todo'):
            path = path[:-4]
        try:
            name = os.path.join(path, root)
            os.mkdir(name)
            return name
        except FileExistsError:
            pass
        try:
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


def remove_if_empty(arg):
    if os.path.isfile(arg) and not os.path.getsize(arg):
        os.remove(arg)


#=======================  RUNNING THE SIMULATION  ==============================

def start_sim(exe, path, args):
    """
    Start simulation in specified working directory (path) as a subprocess
    """
    exe = exe.split(' ')
    outname = os.path.join(path, 'out.txt')
    errname = os.path.join(path, 'err.txt')
    outfile = open(outname, 'w')
    errfile = open(errname, 'w')
    # Option 'cwd=path' sets the current working directory of the subprocess
    if args:
        sub = Popen(exe+args, stdout=outfile, stderr=errfile, cwd=path)
    else:
        sub = Popen(exe, stdout=outfile, stderr=errfile, cwd=path)
    #print(f'Started {exe} process {sub.pid} in `{path}`')
    return sub


def after_sim(sub, path):
    """
    wait for subprocess running in 'path' to complete and return its output
    """
    ret = sub.wait()
    if ret:
        sys.stderr.write(f'{sub.pid} exitcode {ret}')
    else:
        # remove output files if empty:
        remove_if_empty(os.path.join(path, 'out.txt'))
        remove_if_empty(os.path.join(path, 'err.txt'))
    return ret


def info_log(filename, exe, conf, args, pid):
    with open(filename, "w") as f:
        f.write("host      %s\n" % os.getenv('HOSTNAME', 'unknown'))
        f.write("user      %s\n" % os.getenv('USER', 'unknown'))
        f.write("wdir      %s\n" % os.getcwd())
        f.write("exec      %s\n" % exe)
        f.write("args      %s\n" % args)
        f.write("conf      %s\n" % conf)
        f.write("pyid      %s\n" % pid)

def info_start(filename, pid):
    import time
    with open(filename, "a") as f:
        f.write("start     %s\n" % time.asctime())
        f.write("pid       %s\n" % pid)

def info_end(filename, val):
    import time
    with open(filename, "a") as f:
        f.write("status    %s\n" % val)
        f.write("stop      %s\n" % time.asctime())


def run(exe, conf, args, job_name, config_file):
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
    log = os.path.join(wdir, 'log.txt')
    shutil.copyfile(conf, os.path.join(wdir, config_file))
    info_log(log, exe, conf, args, os.getpid())
    sub = start_sim(exe, wdir, args)
    info_start(log, sub.pid)
    ret = after_sim(sub, wdir)
    info_end(log, ret)
    os.chdir(cdir)
    return (ret, wdir)


def start(exe, conf, args, root, config_file):
    """
    Start simulation in a new sub directory, and return immediately.
    The config file `conf` is copied to the sub-directory.
    """
    cdir = os.getcwd()
    if not os.path.isfile(conf):
        raise Error("missing/unreadable config file")
    conf = os.path.abspath(conf)
    wdir = make_directory(root)
    os.chmod(wdir, 504)
    log = os.path.join(wdir, 'log.txt')
    shutil.copyfile(conf, os.path.join(wdir, config_file))
    info_log(log, exe, conf, args, os.getpid())
    sub = start_sim(exe, wdir, args)
    info_start(log, sub.pid)
    return (sub.pid, wdir)


