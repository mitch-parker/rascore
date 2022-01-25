from src.download.modules import *


def compress_output_files(full_path_to_the_file, gzip_mode="on"):
    if gzip_mode == "on":
        with open(Path(full_path_to_the_file), 'rb') as f_in:
            with gzip.open(Path(str(full_path_to_the_file) + '.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
