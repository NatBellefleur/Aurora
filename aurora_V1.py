#ANALYTICAL APPROACH

#phi is the anlge of inclination in degrees of FINCH and therefore also the angle of solar illumination on its panels
#phi is zero at the Northern point of orbit, pi by two at the LTDN and pi at the Southern point.

import math
import scipy.integrate as integrate


def energy_polar(LTDN):
    
    '''four digit analogue time (int) -> energy (float) (watts)
        
        Calculates amount of energy from sun in one revolution given a specific LTDN
    
    ASSUMPTIONS:
    - cycles start at north pole and goes into sunlight direction
    - 100% solar flux efficiency
    - constant orbital velocity
    - no camera movement
    - PERFECTLY POLAR ORBIT'''
    
    #FINCH specific parameters:
    area=0.096576 #square meters
    intensity = 1361 #watts per square meter
    orbit_period=5700 #seconds 
    radius=6371000 #meters
    altitude=550000 #meters
    inclination=97.4*(math.pi/180) #radians 
        
    theta=((LTDN-6)*180/12)*(math.pi/180)
        
    function=lambda phi: math.sin(phi)  #describes changing phi angle
    a=integrate.quad(function,0,math.pi)
    (area_term1,error)=a
        
    max_e=intensity*orbit_period/2*area*math.sin(theta)*area_term1 #max amount of energy in a full cycle (joules)
        
    return max_e

def energy_nonpolar_2(LTDN): 
    
    '''four digit analogue time (int) -> energy (float) (watts)
        
        Calculates amount of energy from sun in one revolution given a specific LTDN
    
    ASSUMPTIONS:
    - cycles start at north pole and goes into sunlight direction
    - 100% solar flux efficiency
    - constant orbital velocity
    - no camera movement'''

    #FINCH specific parameters:
    area=0.096576 #square meters
    intensity= 1361 #watts per square meter
    orbit_period=5700 #seconds 
    radius=6371000 #meters
    altitude=550000 #meters
    inclination=97.4*(math.pi/180) #radians     

    theta=((LTDN-6)*180/12)*(math.pi/180)
    epsilon=abs(math.cos(inclination))*(radius+altitude)
    delta_theta=math.asin(epsilon/radius) #describes the change in epsilon in terms of an angle change
    lower_bound=theta-delta_theta
    upper_bound=theta+delta_theta
        
    function=lambda phi: math.sin(phi) #describes changing phi angle
    a=integrate.quad(function,0,math.pi)
    (area_term1,error)=a

    function=lambda theta: area*math.cos(theta) #describes changing theta angle
    a=integrate.quad(function,lower_bound,upper_bound)
    (area_term2,error)=a

    if area_term1>0 and area_term2>0:
        max_e=intensity*orbit_period/2*area*area_term1*area_term2 #max amount of energy in a full cycle
        
    elif area_term1==0:
        max_e=intensity*orbit_period/2*area*area_term2
        
    elif area_term2==0:
        max_e=intensity*orbit_period/2*area*area_term1
        
    return max_e

#energy_nonpolar_2(12) --> 1.2262055634236836e-12