# MPI converter
# converts a file from, or to MPI version

import sys

if __name__ == "__main__":
    filename = sys.argv[1]
    print(filename)
    f = open(filename)
    print(f)
    txt = f.read()
    txt1 = txt.replace("// MPI COMMENT", "// MPI COMMENT\n  /*")
    txt2 = txt1.replace("// MPI UNCOMMENT", "*/\n  // MPI UNCOMMENT")
    filename_new = filename.replace("_MPI.cpp", ".cpp")
    g = open(filename_new, "w")
    g.write(txt2)
