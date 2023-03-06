import bin
import hapke_photometry
import multiprocessing

def main(obs_dir, rad_dir, band):
    binned = bin.main(obs_dir, rad_dir, band)
    hapke_photometry.main(binned, band)

# def main():
#     if __name__== "__main__":
#         p = multiprocessing.pool()
#         p.map()