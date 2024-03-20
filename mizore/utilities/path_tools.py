import os
def get_subdirectory_paths(root_path):
    paths=[]
    for name in os.listdir(root_path):
        sub_name=os.path.join(root_path,name)
        if os.path.isdir(sub_name):
            paths.append(sub_name)
    return paths