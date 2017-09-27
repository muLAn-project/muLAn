# -*-coding:Utf-8 -*
import sys
sys.path.insert(0, '/Users/Herbie/Recherche/Projets/Miiriads/Softwares/muLAn-project/muLAn/muLAn/models')
import BLcontLLD as BLcontLLD
import fquadBL as fquadBL
import grid_dmcmc as grid_dmcmc
def main():
	models_files = [BLcontLLD]
	models_files.append(fquadBL)
	minim_files = [grid_dmcmc]
	return models_files, minim_files
