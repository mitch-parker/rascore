from src.download.modules import *
from src.download.downloadwithThreadPool import download_with_pool, url_formation_for_pool
from src.renum.shared.SIFTS_tree_parser import SIFTS_tree_parser


def try_SIFTS_tree_parser(default_input_path_to_SIFTS, SIFTS_name):
    product_tree_SIFTS = 0
    for _ in range(3):
        try:
            handle_SIFTS = gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name), 'rt')
            product_tree_SIFTS = SIFTS_tree_parser(handle_SIFTS)
            break
        except EOFError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except ValueError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except OSError:
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except:
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
    return product_tree_SIFTS
