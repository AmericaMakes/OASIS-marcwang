import pandas as pd 
import lxml

from typing import List
import oasis.oasis_setting as oa_s

class PrintFileWriter():

    def __init__(self, filename):
        # TODO init header 
        pass
    
    def write_header(self, h : oa_s.AlsamHeader):
        pass
    
    def write_velocity_profile(self, vl : List[oa_s.AlsamVelocityProfile]):
        pass

    def write_segment_style(self, sl : List[oa_s.AlsamSegmentStyle]):
        pass

    def write_trajectory(self, traj : oa_s.AlsamTrajectory):
        pass
