import os
def clean_examples():

    root_dir = "/Users/khan/MATLAB-Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/"

    for subdir, dirs, files in os.walk(root_dir):
        first_parent = os.path.basename(os.path.normpath(subdir))
        for file in files:
            f_path = os.path.join(subdir, file)
            if ("nodeInfo" in file) | ("MATLABDriveTag" in file):
                os.remove(f_path)
                continue
            if file.endswith('.txt') | file.endswith('.dat') | file.startswith("."):
                print(f_path+"\n")
            else:
                if os.path.isfile(f_path):
                    os.remove(f_path)

        del_empty_dir(subdir)
        # del_empty_dir(dirs)


def clean_conflicts():

    root_dir = "/Users/khan/MATLAB-Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/"

    for subdir, dirs, files in os.walk(root_dir):
        first_parent = os.path.basename(os.path.normpath(subdir))
        for file in files:
            f_path = os.path.join(subdir, file)
            if "conflict_copy" in file:
                os.remove(f_path)

        del_empty_dir(subdir)
        # del_empty_dir(dirs)


def del_empty_dir(dirs):
    if len(os.listdir(dirs)) == 0:
        os.rmdir(dirs)


if __name__ == "__main__":
    clean_conflicts()