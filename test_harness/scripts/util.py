#!/usr/bin/env python

# General useful utilities.



# Simple yes/no question prompted to command line.
# Continues to ask until response in correct format.
# Returns True for yes, False for no
def yes_or_no(question: str):
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} [y,n]: ").strip().lower()
    return (reply == "y")

# Replace all occurences of string with string in file.
def replaceInFile(file,x,y,fileout=None):
    """ replace occurences of x with y in a file """
    if fileout is None: fileout=file
    fd=open(file,'r')
    text=fd.read()
    fd.close()
    fd=open(fileout, 'w')
    text=text.replace(x,y)
    fd.write(text)

# Replace all lines containing 'search' with 'replace' in place
def replaceLineInFile(file,search,replace):
    for line in fileinput.input(file, inplace=True):
        if search in line:
            print(replace)
        else:
            print(line.strip())

# Append line of text to file
def appendToFile(file,text):
    with open(file,"a") as file_object:
        file_object.write(text)

# Replace distutils.util.strtobool
def strtobool(value):
    value = value.lower()
    if value in {"y", "yes", "t", "true", "on", "1"}:
        return True
    if value in {"n", "no", "f", "false", "off", "0"}:
        return False
    raise ValueError(f"invalid truth value {value!r}")
