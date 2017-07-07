#!/usr/bin/python

import numpy as np
import tramp

# CONSTANTS
# in ASTROID coordinates is the X magnet fast (1400 V/s) and the Y one slow
# (57V/s). This leads to moving speed of 31.5 and 1.56 m/s for highest
# energy/momentum (T=225, P=688)
S2US = float(1000000)   # from seconds to microseconds
MS2US = float(1000)     # from milliseconds to microseconds
US2S = float(0.000001)  # from microseconds to seconds
MM2M = float(0.001)     # from mm to m
MP = float(938.27)      # proton rest mass in MeV/c^2

# USER PARAMETERS
maxBeamCurrent = float(2)           # nA
switchEnergyTime = float(2.5*S2US)  # us
signalOnDelay = float(0.5*MS2US)    # us
signalOffDelay = float(0.5*MS2US)   # us
spotSettlingTimeX = float(8*MS2US)  # us
spotSettlingTimeY = float(5*MS2US)  # us
maxMagnetSpeedX = float(3)          # m/s
maxMagnetSpeedY = float(30)         # m/s

scanningMethod = 'step-shoot'    # step-shoot or raster!! (dynamic in the future?)
# print("maxBeamCurrent    = {0}".format(maxBeamCurrent))
# print("switchEnergyTime  = {0}".format(switchEnergyTime))
# print("spotSettlingTimeX = {0}".format(spotSettlingTimeX))
# print("spotSettlingTimeY = {0}".format(spotSettlingTimeY))
# print("maxMagnetSpeedX   = {0}".format(maxMagnetSpeedX))
# print("maxMagnetSpeedY   = {0}".format(maxMagnetSpeedY))
slowDirection = 'x' if maxMagnetSpeedX < maxMagnetSpeedY else 'y'
# print("slowDirection = {0}".format(slowDirection))


def get_magnet_speed_x(momentum):
    """
    Returns magnet speed with proton momentum.
    The speed is proportional to the momentum such that:
    p * v(p) = pmax * v(pmax)
    Thus, knowing pmax and v(pmax), we obtain any speed.
    """
    pmax = 687.67  # MeV/c. if momentum below, scale speed up linearly
    return maxMagnetSpeedX*pmax/momentum


def get_magnet_speed_y(momentum):
    """
    Returns magnet speed with proton momentum.
    The speed is proportional to the momentum such that:
    p * v(p) = pmax * v(pmax)
    Thus, knowing pmax and v(pmax), we obtain any speed.
    """
    pmax = 687.67  # MeV/c. if momentum below, scale speed up linearly
    return maxMagnetSpeedY*pmax/momentum


def get_momentum(energy):
    """Returns proton momentum given energy in MeV."""
    return np.sqrt(energy * (energy + 2*MP))


def get_delivery_time(ngp, max_beam_current):
    """Returns spot delivery time given number of gp and beam current (nA)."""
    return S2US*(ngp*10**9 * 1.6022*10**-19) / (max_beam_current*10**-9)


def get_travel_time(distance, speed):
    """Returns time needed to go from one position to the other"""
    return S2US*abs(distance)*MM2M / speed


def get_timing(e, x, y, w):
    """
    Assigns time and phase to each spot in tramp file.
    
    The implemented model goes as follows:
        1: The driving device (computer) sends a spot to the gantry (assuming already in correct angles).
           This adds a certain signalDelay (750 us).
        2: The magnets move from (0,0) to the spot (x,y). Both at the same time. The total time moving is the slowest.
           The speed of the magnets depend linearly on the momentum with respect to a maximum speed.
        3: The spot position takes some time to settle (5 us, isotropic).
        4: The output energy is changed to the spot energy (constant time of 1 us).
           At this point, we already have the starting delivery time of the spot.
        5: Deliver the spot with a certain beam current (2 nA).
        6: Wait some time before moving on to next spot (1 s).
        7: Repeat 1-6 moving in the changing direction and switching energy when needed until the end of the tramp file.
    
    Between point 5 and 6, we check if the spot ends in the next (or further) phase.
    If it is the case, it is split into two and the phase and ngps are corrected for all affected spots.
    
    There are remaining questions about this model. However, it is intended to be as simple and general as possible:
        - Are the default values reasonable? I find the signal delay specially long.
        - How does the beam current change with the energy/momentum?
        - Are the x and y magnets operated in parallel (implemented) or sequentially?
        - Is the energy switch also operated sequentially with respect to the magnets?
        - Is there a default energy?
    
    In any case, these are sample results of the total counters for reference.
    The tramp file employed contained 1454.737275 gp in 60516 spots
        Number of layers   =    36
        Number of spots    = 60516
        Number of gps      =  1454.734
        Moving beam        =   315.295 s
        Switching Energies =    90.000 s
        Delivery time      =   116.539 s
        Signal Delay time  =    60.516 s
        --------------------------------
        Total treatment    =   582.350 s
    
    In a more realistic case:
        Number of layers   =   17
        Number of spots    = 1741
        Number of gps      =  375.671
        Moving beam        =   16.174 s
        Switching Energies =   42.500 s
        Delivery time      =   30.095 s
        Signal Delay time  =    1.679 s
        -------------------------------
        Total treatment    =   90.448 s
    """

    source = tramp.Tramp(e=e, x=x, y=y, w=w)

    # TOTAL COUNTERS
    total_time_moving_beam = float(0)
    total_time_switching_e = float(0)
    total_time_delivering = float(0)
    total_time_signal_delay = float(0)
    
    previous_x = float(0)
    previous_y = float(0)
    previous_e = float(-1)
    time = float(0)

    # Loop over each spot
    i = int(0)
    while i < len(source.spots):
        # If new layer, go to (0,0) and switch energy
        spot = source.spots[i]

        if spot.energy != previous_e:
            p = get_momentum(spot.energy)
            travel_time_x = get_travel_time(previous_x - 0, get_magnet_speed_x(p)) +\
                spotSettlingTimeX * (previous_x != 0)
            travel_time_y = get_travel_time(previous_y - 0, get_magnet_speed_y(p)) +\
                spotSettlingTimeY * (previous_y != 0)
            previous_x = 0
            previous_y = 0
            time += switchEnergyTime + max(travel_time_x, travel_time_y)
            total_time_switching_e += switchEnergyTime
            total_time_moving_beam += max(travel_time_x, travel_time_y)
        
        # Wait signalOnDelay while beam builds-up if not split spot
        time += signalOnDelay
        total_time_signal_delay += signalOnDelay
        
        # move beam to spot in layer
        travel_time_x = 0
        travel_time_y = 0
        if spot.x != previous_x:
            travel_time_x = get_travel_time(spot.x - previous_x, get_magnet_speed_x(p)) + spotSettlingTimeX
        if spot.y != previous_y:
            travel_time_y = get_travel_time(spot.y - previous_y, get_magnet_speed_y(p)) + spotSettlingTimeY
        time += max(travel_time_x, travel_time_y)
        total_time_moving_beam += max(travel_time_x, travel_time_y)
        
        # assign starting values to current spot
        spot.time = float(time)

        # add delivery time if not repeated spot
        spot_delivery_time = get_delivery_time(spot.ngp, maxBeamCurrent)
        time += spot_delivery_time
        total_time_delivering += spot_delivery_time

        # print(spot)
        
        # Add doseDelay if spot is not repeated.
        # We don't care if it's the last in treatment
        time += signalOffDelay
        total_time_signal_delay += signalOffDelay
        
        # Prepare next spot
        previous_x = spot.x
        previous_y = spot.y
        previous_e = spot.energy
        i += 1
    
    # print("===== Summary (s) ==========")
    # print(" Number of layers   = {}".format(source.get_number_layers()))
    # print(" Number of spots    = {}".format(source.number_spots))
    # print(" Moving beam        = {0:.3f}".format(round(US2S*total_time_moving_beam, 3)))
    # print(" Switching Energies = {}".format(US2S*total_time_switching_e))
    # print(" Delivery time      = {0:.3f}".format(round(US2S*total_time_delivering, 3)))
    # print(" Signal Delay time  = {}".format(US2S*total_time_signal_delay))
    # print(" ---------------------------")
    # print(" Total treatment    = {0:.3f}".format(round(US2S*time, 3)))
    # print("============================")

    summary = """ Summary (s):
 Moving       = {:.3f}
 Energies     = {}
 Delivery     = {:.3f}
 Signal Delay = {:.3f}
 ---------------------
 Total        = {:.3f}""".format(
        round(US2S*total_time_moving_beam, 3),
        US2S*total_time_switching_e,
        round(US2S*total_time_delivering, 3),
        US2S*total_time_signal_delay,
        round(US2S*time, 3))

    return np.array([US2S*source.spots[t].time for t in range(len(source.spots))]), summary
