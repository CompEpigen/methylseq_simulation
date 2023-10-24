import numpy as np
import random, logging

def set_seed(seed) -> None:
	random.seed(seed)
	np.random.seed(seed)

def set_logger() -> None:
	logging.basicConfig(format='%(asctime)s:%(name)s:%(levelname)s:\n\t%(message)s')

	