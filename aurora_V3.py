#VECTOR APPROACH --> GENERAL SATELLITE

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from scipy.integrate import simps
from scipy.spatial.transform import Rotation


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

def sunlit(panels:dict, rotations:Rotation, sun_vectors:np.ndarray, timesteps: np.ndarray, ECI_vectors: np.ndarray):

  if np.any(panels.values() == np.array([1, 0, 0])):
    #we are applying rotations to "ECI x" to produce "spacecraft x"
    spacecraft_orientations_xpos = rotations.apply(np.array([1, 0, 0])) #takes 1D slices of the rotations array on the given x axis. 
    alignments_xpos = np.sum(spacecraft_orientations_xpos * sun_vectors, axis=1) #sums the array elements of the spacecraft oriention and the sunlight vectors in the chosen x axis. 

  elif np.any(panels.values() == np.array([-1, 0, 0])):
    #we are applying rotations to "ECI -x" to produce "spacecraft x-"
    spacecraft_orientations_xneg = rotations.apply(np.array([-1, 0, 0])) 
    alignments_xneg = np.sum(spacecraft_orientations_xneg * sun_vectors, axis=1) 

  elif np.any(panels.values() == np.array([0, 1, 0])):
    #we are applying rotations to "ECI y" to produce "spacecraft y"
    spacecraft_orientations_ypos = rotations.apply(np.array([0, 1, 0])) 
    alignments_ypos = np.sum(spacecraft_orientations_ypos * sun_vectors, axis=1)
    
  elif np.any(panels.values() == np.array([0, -1, 0])):
    #we are applying rotations to "ECI -y" to produce "spacecraft -y"
    spacecraft_orientations_yneg = rotations.apply(np.array([0, -1, 0])) 
    alignments_yneg = np.sum(spacecraft_orientations_yneg * sun_vectors, axis=1)

  elif np.any(panels.values() == np.array([0, 0, 1])):
    #we are applying rotations to "ECI z" to produce "spacecraft z"
    spacecraft_orientations_zpos = rotations.apply(np.array([0, 0, 1])) 
    alignments_zpos = np.sum(spacecraft_orientations_zpos * sun_vectors, axis=1)

  elif np.any(panels.values() == np.array([0, 0, -1])):
    #we are applying rotations to "ECI -z" to produce "spacecraft -z"
    spacecraft_orientations_zneg = rotations.apply(np.array([0, 0, -1])) 
    alignments_zneg = np.sum(spacecraft_orientations_zneg * sun_vectors, axis=1)


  zeroes = np.zeros(timesteps.shape[0]) #creates an array containing only zero values

                                        
  alignments_xpos[np.where(alignments_xpos < 0)] = zeroes[np.where(alignments_xpos < 0)]   #replace the negative alignment values with the specific zeros in the mirrored position in the zeros array                                                                            # alignments < 0 == True/False -- Only activated if true
  xpos = shadowed(ECI_vectors,sun_vectors)
  alignments_xpos[np.where(xpos)] = zeroes[np.where(xpos)]                               #replacing the True responses (spacecraft in the shadow) with zero place holders from the zeros array. 


  alignments_xneg[np.where(alignments_xneg < 0)] = zeroes[np.where(alignments_xneg < 0)]                                                                       
  xneg = shadowed(ECI_vectors,sun_vectors)
  alignments_xneg[np.where(xneg)] = zeroes[np.where(xneg)]   


  alignments_ypos[np.where(alignments_ypos < 0)] = zeroes[np.where(alignments_ypos < 0)]                                                                    
  ypos = shadowed(ECI_vectors,sun_vectors)
  alignments_ypos[np.where(ypos)] = zeroes[np.where(ypos)]  


  alignments_yneg[np.where(alignments_yneg < 0)] = zeroes[np.where(alignments_yneg < 0)]                                                                             
  yneg = shadowed(ECI_vectors,sun_vectors)
  alignments_yneg[np.where(yneg)] = zeroes[np.where(yneg)]   


  alignments_zpos[np.where(alignments_zpos < 0)] = zeroes[np.where(alignments_zpos < 0)]                                                                    
  zpos = shadowed(ECI_vectors,sun_vectors)
  alignments_zpos[np.where(zpos)] = zeroes[np.where(zpos)]  


  alignments_zneg[np.where(alignments_zneg < 0)] = zeroes[np.where(alignments_zneg < 0)]                                                                             
  zneg = shadowed(ECI_vectors,sun_vectors)
  alignments_zneg[np.where(zneg)] = zeroes[np.where(zneg)]  

  alignments = alignments_xpos + alignments_xneg + alignments_ypos + alignments_yneg + alignments_zpos + alignments_zneg
  x = panels.keys()
  areas = sum(x)

  results = [alignments, areas]
  return results


def aurora (panels,rotations: Rotation, sun_vectors: np.ndarray, times: np.ndarray, ECI_vectors: np.ndarray, intensity: float)-> float:

  results = sunlit(panels,rotations,sun_vectors,times,ECI_vectors)

  total_energy = simps (intensity * results[1] * results[0],times)

  average_power = total_energy/times[-1]

  return average_power


# COMMAND WINDOW

# finch_panels = {0.024144:np.array([1, 0, 0]), 0.024144:np.array([-1, 0, 0]), 0.024144:np.array([0, 1, 0]), 0.024144:np.array([0, -1, 0])}                                                                                                                                      
# intensity = 1361

# x = aurora (finch_panels,rotations,sun_vectors,times,ECI_vectors,intensity)

# display (x)

#next work on CLASS APPROACH --> GENERAL SATELLITE