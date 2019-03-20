# Blender rendering

The startup files ‘actin.blend’ and ‘cytosim.blend’ are identical.
You should keep one as a backup.

Scripts to import Cytosim's objects info Blender:

    import-actin.py    add monomers from files ‘cymart/actin####.txt’
    import-links.py    add links from files ‘cymart/link####.txt’
    delete-actin.py    delete stuff added by the other two scripts

On cytosim’s side, the tool ‘cymart’ correctly makes these files.
Check the help, but essentially, it should be called like this:

    cymart actin style=actin

The output files should then be copied in a subfolder ‘cymart’ of your blender directory.

# Python error messages

On MacOS, Python error messages are sent to the terminal, if you start blender from a terminal. So open a Terminal window, and enter somthing like this:

    /Applications/blender-2.79/blender.app/Contents/MacOS/blender actin.blend

You will then be able to see the error messages.

# More

There is another version of import-actin.py:

    import-actin-ico.py    Version that uses a icosahedral mesh
