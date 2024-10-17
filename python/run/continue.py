#!/usr/bin/env python3
#
# Continue a batch of unfinished simulations
#
# Helen Saville, Cambridge University 2024
"""

------- Description ------
Continues simulations not finished on the cluster, until the final time specified in their original config file

'run' stage:
Prepares unfinished simulations to be continued:
 - creates new folders with the same naming pattern r%%%% as the old folders
 - extracts final frame of unfinished simulation and moves it to the new folder,
 - copies properties.cmi to new folder,
 - creates config file to continue the simulation.
 Continues unfinished simulations:
 - either using submit_one.py if on the cluster, or scan.py if on local machine

 'combine' stage:
 Combines original and continued data:
 - everything copied to the new folder
 - objects.cmo files concatenated (but objects_original.cmo and objects_cont.cmo also retained)

 'clean' stage:
 Deletes folders and files to do with unfinished simulations, retaining only combined results
 Everything ends up in ORIGINAL_RUN_DIRECTORY_cont/save_cont/r%%%%

'check' stage:
Check if continued simulations have run to the end of the original run time
If they have not, move them into ORIGINAL_RUN_DIRECTORY_cont/r%%%%, so that can repeat the process again

 Run this script in a directory that contains
 1) the original run directory, and 2) the executables 'sim', 'frametool', and 'report'
Also, the directory that continue.py is in should also contain submit_one.py and scan.py.

Needs python version >=3.7

------- Syntax ---------
continue.py original_run_directory [exe=EXECUTABLE] [stage=STR] [cluster=BOOL] [jobs=INT]

compulsory:
original_run_directory: directory containing unfinished simulations, whose folder names start with r%%%%
stage: run, combine, clean, or check
exe: path to executable, relative to current folder

optional:
cluster = true/false (default=true): Run on cluster (SLURM) or local machine
jobs = INT (default=1): if running on local machine, set number of jobs


"""

import sys, shutil, os, re, subprocess


def make_new_config(frames_run, time_run, time_step, run_steps, nb_frames):
    # make a new config

    run_left = int(run_steps - time_run/time_step)
    frames_left = nb_frames - frames_run
    config = open('config.cym', 'w')
    text = "read properties.cmp \nimport all objects.cmi  \nrun %s system {  nb_frames = %s }"
    config.write(text % (run_left, frames_left))
    frames_txt = open('frames_left.txt', 'w')
    frames_txt.write(str(frames_left))


def get_final_frame(frames_run, frametool):
    # extract final frame to objects.cmi

    final_frame = frametool + ' objects.cmo ' + str(frames_run) + ' > objects.cmi'
    os.system(final_frame)


def process_old_config():
    # extract the time step, number of run steps, and number of frames from the original config file
    #     (assumes there is space around equal signs)

    config = open('config.cym', 'r')
    time_step = 0
    run_steps = 0
    nb_frames = 0
    for line in config:
        s = line.split()
        if len(s) > 1:
            if re.match('time_step', s[0]):
                time_step = s[2]
            elif re.match('run', s[0]):
                run_steps += int(s[1])
            elif re.match('nb_frames', s[0]):
                nb_frames = s[2]
    return float(time_step), int(run_steps), int(nb_frames)


def is_float(string):
    if string.replace(".", "").isnumeric():
        return True
    else:
        return False


def get_run_data(frametool):
    # get number of frames and run time completed during the original run
    time_run = 0

    command = [frametool, 'last', 'objects.cmo']
    time_run_data = subprocess.run(command, capture_output=True, text=True, check=True, encoding='latin-1')
    time_run_data = time_run_data.stdout
    lines = time_run_data.splitlines()
    for line in lines:
        s = line.split()
        if len(s) > 1:
            if re.match('#time', s[0]):
                time_run = float(s[1])
                break

    command = [frametool, 'objects.cmo']
    nb_frames = subprocess.run(command, capture_output=True, text=True, check=True)
    nb_frames = nb_frames.stdout
    nb_frames = nb_frames.strip()
    nb_frames = nb_frames.split()
    frames_run = int(nb_frames[1]) - 1  # frametool counts the initialisation frame in the number of frames

    return frames_run, time_run


def prepare_one(run_path, run_cont_path, frametool):
    # prepare new folders to continue simulations in

    os.chdir(run_path)
    # if steps > 0:
    frames_run, time_run = get_run_data(frametool)
    time_step, run_steps, nb_frames = process_old_config()
    if not os.path.isfile('objects.cmi'):
        get_final_frame(frames_run, frametool)

    shutil.move(run_path+'/objects.cmi', run_cont_path+'/objects.cmi')
    shutil.copy(run_path+'/properties.cmp', run_cont_path+'/properties.cmp')

    os.chdir(run_cont_path)
    make_new_config(frames_run, time_run, time_step, run_steps, nb_frames)


def prepare(run_paths, run_cont_paths, frametool):
    for run_path, run_cont_path in zip(run_paths, run_cont_paths):
        if not os.path.isdir(run_cont_path):
            os.mkdir(run_cont_path)
        else:
            print('a continuation directory already exists, so continued simulations may have already been run')
            return 0
        # processed unfinished simulation data and prepare to continue
        prepare_one(run_path, run_cont_path, frametool)
    return 1


def run(run_cont_paths, cluster, jobs, exe, runmod):
    # run unfinished simulations in the new folders
    args = [exe]
    args.extend(run_cont_paths)
    if cluster:
        runmod.main(args)
    else:
        args.append('jobs='+jobs)
        runmod.main(args)


def combine(run_paths, run_cont_paths, run_names, save_paths, frametool):
    # combine data for original and continued runs, used in main()

    # loop through original and continued simulations, and move everything to continued simulation folder
    for (run_name, run_path, run_cont_path, save_path) in zip(run_names, run_paths, run_cont_paths, save_paths):

        # check if results for this folder have already been combined and moved
        if not os.path.isdir(save_path):
            # check that continuation simulation created objects.cmo
            if os.path.isfile(os.path.join(run_cont_path, 'objects.cmo')):
                # rename config and objects in continued simulation folder
                os.rename(run_cont_path+'/config.cym', run_cont_path+'/config_cont.cym')
                os.rename(run_cont_path+'/objects.cmo', run_cont_path+'/objects_cont.cmo')

                # copy everything in original simulation folder to continued simulation folder
                for f in os.listdir(run_path):
                    src = os.path.join(run_path, f)
                    dst = os.path.join(run_cont_path, f)
                    shutil.copy(src, dst)

                # rename original objects.cmo in continued simulation folder
                try:
                    os.rename(run_cont_path+'/objects.cmo', run_cont_path+'/objects_orig.cmo')
                except:
                    print('no original objects.cmo in %s' % run_path)

                # remove first frame from objects_cont.cmo, as it is the same as the final frame in objects_orig.cmo
                os.chdir(run_cont_path)

                read_frames_left = (open('frames_left.txt', 'r')).read()
                frames_left = int(read_frames_left)
                frames = frametool + ' objects_cont.cmo ' + '1:' + str(frames_left) + ' > objects_cont1.cmo'
                os.system(frames)

                # concatenate results from original and continued simulations
                os.system('cat objects_orig.cmo objects_cont1.cmo > objects.cmo')

                # move everything from path_cont/r%%%% to path_cont/save_cont/r%%%%
                shutil.move(run_cont_path, save_path)

            else:
                print('No objects.cmo in continuation of %s' % run_name)
        else:
            print('Results for %s have already been combined and moved to save_cont' % run_name)


def clean(run_paths, save_paths):
    for run_path, save_path in zip(run_paths, save_paths):
        # check if results for this folder have already been combined and moved
        if os.path.isdir(save_path):
            # remove all files and folders associated with unfinished simulation
            shutil.rmtree(run_path)
            os.chdir(save_path)
            os.remove('objects_orig.cmo')
            os.remove('objects_cont.cmo')
            os.remove('objects.cmi')
            os.remove('config_cont.cym')
            os.remove('objects_cont1.cmo')
            os.remove('frames_left.txt')


def check(run_names, save_paths, frametool):
    for run_name, save_path in zip(run_names, save_paths):
        os.chdir(save_path)
        # extract actual simulation time
        time_step, run_steps, nb_frames = process_old_config()
        # extract desired simulation time
        frames_run, time_run = get_run_data(frametool)
        frames_left = int(nb_frames - frames_run)
        if frames_left != 0:
            shutil.move(save_path, '../../')
            print('need to continue %s for %s frames' % (run_name, frames_left))


def list_paths(path):
    # get list of directories for original simulations, continued simulations, and combined results

    path_cont = path + "_cont"

    # make directory for files with which to continue simulations
    if not os.path.isdir(path_cont):
        os.mkdir(path_cont)

    # absolute paths
    cwd = os.getcwd()
    path = os.path.join(cwd, path)
    path_cont = os.path.join(cwd, path_cont)
    save_path = os.path.join(path_cont, 'save_cont')

    # pattern for directories that start with 'r' followed by an integer
    pattern = re.compile(r'^r\d+')

    # get path lists of unfinished simulation run folder and folders to continue them in
    run_paths = []
    run_cont_paths = []
    run_names = []
    save_paths = []
    for subdir in os.listdir(path):
        # check if subdirectory has the naming pattern r%%%% of a simulation run
        matched = pattern.match(subdir)
        if matched:
            name = matched.group(0)
            run_names.append(name)
            run_paths.append(os.path.join(path, subdir))
            run_cont_paths.append(os.path.join(path_cont, name))
            save_paths.append(os.path.join(save_path, name))
    return save_path, run_paths, run_cont_paths, run_names, save_paths


def check_true(arg):
    if arg == 'False' or arg == 'false' or arg == '0':
        return False
    else:
        return True


def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)


def main(args):

    path = ''
    cluster = True
    jobs = '1'
    stage = ''
    exe = os.path.abspath('sim')

    # parse arguments:
    for arg in args:
        if os.path.isdir(arg):
            path = arg
        elif executable(arg):
            exe = os.path.abspath(arg)
        else:
            [key, equal, val] = arg.partition('=')
            if key == 'path':
                path = val
            elif key == 'stage':
                stage = val
            elif key == 'cluster':
                cluster = check_true(val)
            elif key == 'jobs':
                jobs = val
            elif key == 'exe':
                exe = os.path.abspath(val)

    if cluster:
        import scan
        runmod = scan
    else:
        import submit_one
        runmod = submit_one

    save_path, run_paths, run_cont_paths, run_names, save_paths = list_paths(path)

    frametool = os.path.abspath('frametool')


    # prepare and run continuation sims
    if stage == 'run':

        # Loop through subdirectories to find unfinished runs prepare to continue them
        do_run = prepare(run_paths, run_cont_paths, frametool)

        # Continue each unfinished simulation
        if do_run:
            run(run_cont_paths, cluster, jobs, exe, runmod)

    # Combine original and continued data
    elif stage == 'combine':
        if not os.path.isdir(save_path):
            os.mkdir(save_path)
        combine(run_paths, run_cont_paths, run_names, save_paths, frametool)

    # Cleanup unneeded files
    elif stage == 'clean':
        if os.path.exists(save_path) and os.listdir(save_path):
            clean(run_paths, save_paths)
        else:
            print('no results to cleanup yet: need to run and/or combine first')

    # Check if continued simulations ran to the full simulation time
    elif stage == 'check':
        if os.path.exists(save_path) and os.listdir(save_path):
            check(run_names, save_paths, frametool)
        else:
            print('no results to check yet: need to run and/or combine first')

    else:
        print('stage not specified')


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "You must specify a directory containing simulations to continue (for instructions, "
            "invoke with option '--help')")
    elif sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])
