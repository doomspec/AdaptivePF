import pathlib
import json
import pickle
import os


def put_log_item_in_json(path, keyname, value):
    pathlib.Path(path).touch()
    with open(path, "r") as f:
        raw = f.read()
    if len(raw) != 0:
        log_dict = json.loads(raw)
    else:
        log_dict = dict()
    log_dict[keyname] = value
    with open(path, "w") as f:
        json.dump(log_dict, f)


def get_paths_by_heading(folder_path, heading):
    """
    example: H2O_xxxxx -> heading = H2O
    """
    paths = []
    g = os.walk(folder_path)
    for path, dir_list, _file_list in g:
        for dir_name in dir_list:
            sub_path = os.path.join(path, dir_name)
            file_name = dir_name.split("_")
            if file_name[0] == heading:
                paths.append(sub_path)
    return paths
