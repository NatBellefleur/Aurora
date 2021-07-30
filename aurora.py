import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from scipy.integrate import simps
from scipy.spatial.transform import Rotation

#Testing push to github

def shadowed(sat_ECI: np.ndarray, sun_vec: np.ndarray, R_Earth = 6378, R_Sun = 149598023)-> bool:

    mag_sat_ECI = np.linalg.norm(sat_ECI, axis=-1)      
   
    theta = np.arccos(np.sum(sat_ECI * sun_vec, axis=-1)/mag_sat_ECI)
      # Angle between Sat and Sun

    theta_1 = np.arccos(R_Earth / mag_sat_ECI)
      # Angle between tangent and sun

    theta_2 = np.arccos(R_Earth / R_Sun)
      # Angle between tangent and spacecraft

    return theta_1 + theta_2 < theta
    #True - satellite in Earth's shadow
    #False - satellite in sunlight

def aurora (rotations: Rotation, sun_vectors: np.ndarray, timesteps: np.ndarray, ECI_vectors: np.ndarray, area: float, intensity = 1361.0)-> float:

    spacecraft_orientations = rotations.apply(np.array([0, 0, 1]))
    alignments = np.sum(spacecraft_orientations * sun_vectors, axis=1)

    shadowed(ECI_vectors, sun_vectors)

    zeroes = np.zeros(timesteps.shape[0])
    alignments[np.where(alignments < 0)] = zeroes[np.where(alignments < 0)]

    x = shadowed(ECI_vectors,sun_vectors)
    alignments[np.where(x)] = zeroes[np.where(x)]

    total_energy = simps (intensity * area * alignments,timesteps)

    average_power = total_energy/timesteps[-1]

    return average_power

def plot_power (timesteps, alignments):

    plt.figure(figsize=(12,7))
    plt.plot(timesteps[:2000], alignments[:2000])
    plt.show()