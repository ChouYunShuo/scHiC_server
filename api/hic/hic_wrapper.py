
def parse_hic_path(s: str):
    # /path/to/file::/path/to/grp

    parts = s.split("::")

    file_path, group_path = parts
    if not group_path.startswith("/"):
        group_path = "/" + group_path

    return file_path, group_path


class Hic_wrapper(object):
    def __init__(self, hic_path: str):
        self.filepath, self.rootpath = parse_hic_path(hic_path)


"""
functionalities

1. given a .hic file, generate a single cell multi-resolution h5 file

2. data query 

3. API for frontend










"""
