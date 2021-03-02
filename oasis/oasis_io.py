import pandas as pd 
import lxml

from typing import List, Dict
import oasis.oasis_setting as oa_s
from datetime import date

class PrintFileLoader():
    def __init__(self, filename):
        self.filename = filename

    def load_header(self) -> oa_s.AlsamHeader:
        df = pd.read_excel(self.filename, sheet_name= 1, engine="xlrd")
        layer_thickness = df.iloc[4, 2]
        return oa_s.AlsamHeader(date.today(), layer_thickness)
    
    def load_velocity_profile(self) -> Dict[str, oa_s.AlsamVelocityProfile]:
        df = pd.read_excel(self.filename, 
                sheet_name=2, index_col= 0, header = 5, usecols="A:H")
        vel_alsam = {}
        for index, row in df.iterrows():
            vel_alsam[index] = oa_s.AlsamVelocityProfile(index, row.iloc[0], 
                row.iloc[2], row.iloc[3], row.iloc[4], row.iloc[5], row.iloc[6], row.iloc[1])
        return vel_alsam
    
    def load_segment_style(self) :
        df = pd.read_excel(self.filename, sheet_name=3, index_col= 0, header = 7, usecols="A:J")
        seg_style = {}
        for index, row in df.iterrows():
            traveler = oa_s.AlsamTraveler(row.iloc[1], row.iloc[2], row.iloc[3])
            seg_style[index] = oa_s.AlsamSegmentStyle(index, row.iloc[0], traveler)

        return seg_style

    def load_parts(self):
        regions = pd.read_excel(self.filename, sheet_name=4, header=5, usecols="A:O", index_col=0)
        parts = pd.read_excel(self.filename, sheet_name = 5, header=5, usecols="A:G")

        regions_d = {}
        for index , row in regions.iterrows():
            p_r = oa_s.PrintRegions(index,*row.iloc[0:13].values.tolist())
            regions_d[index] = p_r

        part_d = []
        for index, row in parts.iterrows():
            part = oa_s.PartsSetting(*row.iloc[0:7].values.tolist())
            part_d.append(part)
        
        return regions_d, part_d




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
