"""
Reorginize LAC product in clean file structure.
"""
import os
import sys


def __poststeps__():
    if len(sys.argv) != 2:
        print('Usage: postproduce.py OUTPUT_PATH')
        sys.exit(1)
    poststeps(sys.argv[1])


def poststeps(base_folder):
    """simple cleanup steps."""
    tmp_path = os.path.join(base_folder, 'TMP')
    if not os.path.exists(tmp_path):
        print('Output folder does not exists')
        sys.exit(1)
    for path in os.listdir(tmp_path):
        src_path = os.path.join(tmp_path, path)
        if path.endswith('.TIF') or path.endswith('.jpg'):
            os.rename(src_path, os.path.join(base_folder, 'GRANULE', path))
        else:
            os.rename(src_path, os.path.join(base_folder, 'AUX_DATA', path))
    os.rmdir(tmp_path)


if __name__ == '__main__':
    __poststeps__()
