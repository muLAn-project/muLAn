# -*-coding:Utf-8 -*
# -----------------------------------------------------------------------------
# This file will be removed in a next version - kept for compatibility
# -----------------------------------------------------------------------------
try:
    import sys
    import os
except:
    print "ERROR: The standard packages 'sys' and 'os' are required."

#   Initialisation without any packages loaded but sys and os
# -----------------------------------------------------------------------------
full_path = os.path.realpath(__file__)
text = full_path.split('/')
a = ''
i = 0
while i < len(text)-1:
   a = a + text[i] + '/'
   i = i + 1
full_path = a

full_path_code = ''
path_lib_ext = 'None'

filename = full_path + 'advancedsetup.ini'
file = open(filename, 'r')
for line in file:
    line = line.replace('\n', '')
    line = [a.strip() for a in line.split(':')]

    if len(line) > 1:
        if ((line[0])[0] != '#') & ((line[0])[0] != ';')  & (line[1] != ''):

            # External packages
            if line[0] == 'PathLocalPythonPackages':
                path_lib_ext = line[1]
                if path_lib_ext[-1]!='/': path_lib_ext = path_lib_ext + '/'
                sys.path.insert(0, path_lib_ext)

            # Path to the code
            if line[0] == 'Code':
                full_path_code = line[1]
                if full_path_code[-1]!='/': full_path_code = full_path_code + '/'
file.close()

if full_path_code == '': sys.exit('I do not find a path for muLAn.')

#   Load the logo
# -----------------------------------------------------------------------------
logo = ''
filename = full_path_code + 'Logo.txt'
file = open(filename, 'r')
for line in file:
    logo = logo + line
file.close()

logo = logo + '\033[2mMicro-Lensing Analysis software\033[0m\n'
logo = logo + '\033[2mVersion '

filename = full_path_code + 'Version.txt'
file = open(filename, 'r')
for line in file:
    logo = logo + line
    break
file.close()

logo = logo + '\033[0m'

print logo

#   Test of the required packages
# -----------------------------------------------------------------------------
packages_required = []
filename = full_path_code + 'PackageVersions.txt'
file = open(filename, 'r')
for line in file:
    line = line.replace('\n', '')
    line = [a.strip() for a in line.split(',')]

    if len(line) > 1:
        packages_required.append(line)
file.close()

text = 'import sys\n'
for i in xrange(len(packages_required)):
    if (packages_required[i][1] == 'None') | (packages_required[i][1] == ''):
        text = text\
                + 'try:\n'\
                + '    import ' + packages_required[i][0] + '\n'\
                + 'except:\n'\
                + "    sys.exit('The required package " + packages_required[i][0]\
                + " is not installed.')\n"
    else:
        text = text\
                + 'try:\n'\
                + '    import ' + packages_required[i][0] + '\n'\
                + '    version = ' + packages_required[i][0] + '.__version__\n'\
                + '    version2 = [int(a) for a in version.split(".")]\n'\
                + '    version_ref = "' + packages_required[i][1] + '"\n'\
                + '    version_ref2 = [int(a) for a in version_ref.split(".")]\n'\
                + '    text_retard = "User Warning: the version of '\
                + packages_required[i][0] + ' should be >= '\
                + packages_required[i][1] + '. The active version is " + version + "."\n'\
                + '    text_avance = "User Warning: the optimal version of '\
                + packages_required[i][0] + ' is '\
                + packages_required[i][1] + '. The active version is " + version + "."\n'\
                + '    if version_ref2[0]>version2[0]: print text_retard\n'\
                + '    elif version_ref2[0]<version2[0]: print text_avance\n'\
                + '    else:\n'\
                + '        if version_ref2[1]>version2[1]: print text_retard\n'\
                + '        elif version_ref2[1]<version2[1]: print text_avance\n'\
                + '        else:\n'\
                + '            if version_ref2[2]>version2[2]: print text_retard\n'\
                + '            elif version_ref2[2]<version2[2]: print text_avance\n'\
                + 'except:\n'\
                + "    sys.exit('The required package " + packages_required[i][0]\
                + " is not installed.')\n"

filename = full_path + '.packages_check.py'
file = open(filename, 'w')
file.write(text)
file.close()
execfile(filename)
os.remove(filename)
# ----------------------------------------------------------------------
#   Packages
# ----------------------------------------------------------------------
import subprocess
# ----------------------------------------------------------------------
#   Functions
# ----------------------------------------------------------------------
def bash_command(text):
    proc = subprocess.Popen(text, shell=True, executable="/bin/bash")
    proc.wait()
# ----------------------------------------------------------------------
def start(full_path, full_path_code, path_lib_ext):
    # Python local libraries
    filename = full_path + '.pythonexternallibpath'
    file = open(filename, 'w')
    file.write(path_lib_ext)
    file.close()

    filename = full_path_code + '.pythonexternallibpath'
    file = open(filename, 'w')
    file.write(path_lib_ext)
    file.close()

    # Run code
    filename = full_path_code + 'main.py'
    text = 'python ' + filename + ' ' + full_path
    bash_command(text)

    # Remove files
    filename = full_path + '.pythonexternallibpath'
    if os.path.exists(filename): os.remove(filename)
    filename = full_path_code + '.pythonexternallibpath'
    if os.path.exists(filename): os.remove(filename)
    filename = full_path + '.emergencystop'
    if os.path.exists(filename): os.remove(filename)

    # Free memory
    del filename, text, file
# ----------------------------------------------------------------------
def stop(full_path):
    # Safe Emergency Stop
    if os.path.exists(full_path + '.emergencystop'):
        os.remove(full_path + '.emergencystop')
        file = open(full_path + '.emergencystop', 'w')
        file.write('1')
        file.close()
# ----------------------------------------------------------------------
def order(full_path, full_path_code, path_lib_ext):
    # Python local libraries
    filename = full_path + '.pythonexternallibpath'
    file = open(filename, 'w')
    file.write(path_lib_ext)
    file.close()

    filename = full_path_code + '.pythonexternallibpath'
    file = open(filename, 'w')
    file.write(path_lib_ext)
    file.close()

    # Run stand-alone package order_ChainsResults
    filename = full_path_code + 'packages/order_ChainsResults.py'
    text = 'python ' + filename + ' ' + full_path
    bash_command(text)

    # Remove files
    # filename = full_path + '.pythonexternallibpath'
    # if os.path.exists(filename): os.remove(filename)
    # filename = full_path_code + '.pythonexternallibpath'
    # if os.path.exists(filename): os.remove(filename)
    # filename = full_path + '.emergencystop'
    # if os.path.exists(filename): os.remove(filename)

    # Free memory
    del filename, text, file
# ----------------------------------------------------------------------
if (__name__ == "__main__"):
    if len(sys.argv) > 1:
        todo = sys.argv[1]
        if todo == 'start':
            start(full_path, full_path_code, path_lib_ext)
        elif todo == 'stop':
            stop(full_path)
        elif todo == 'order':
            order(full_path, full_path_code, path_lib_ext)
    else:
        print "Please use keyword 'start' or 'stop'. 'Start' mode is used."
        start(full_path, full_path_code, path_lib_ext)
