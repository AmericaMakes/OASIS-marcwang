import pandas as pd
from lxml import etree

from typing import List, Dict
import oasis.oasis_setting as oa_s
from datetime import date

from shapely.geometry import Polygon, MultiLineString


class PrintFileLoader():
    def __init__(self, filename):
        self.filename = filename
        self.header = self.load_header()
        self.vel_profile = self.load_velocity_profile()
        self.seg_style = self.load_segment_style()

        self.regions, self.parts = self.load_parts()

    def load_header(self) -> oa_s.AlsamHeader:
        df = pd.read_excel(self.filename, sheet_name=1, engine="xlrd")
        layer_thickness = df.iloc[4, 2]
        return oa_s.AlsamHeader(date.today(), layer_thickness)

    def load_velocity_profile(self) -> Dict[str, oa_s.AlsamVelocityProfile]:
        df = pd.read_excel(self.filename,
                           sheet_name=2, index_col=0, header=5, usecols="A:H")
        vel_alsam = {}
        for index, row in df.iterrows():
            vel_alsam[index] = oa_s.AlsamVelocityProfile(index, row.iloc[0],
                                                         row.iloc[2], row.iloc[3], row.iloc[4], row.iloc[5], row.iloc[6], row.iloc[1])
        return vel_alsam

    def load_segment_style(self) -> Dict[str, oa_s.AlsamSegmentStyle]:
        df = pd.read_excel(self.filename, sheet_name=3,
                           index_col=0, header=7, usecols="A:J")
        seg_style = {}
        for index, row in df.iterrows():
            traveler = oa_s.AlsamTraveler(
                row.iloc[1], row.iloc[2], row.iloc[3])
            seg_style[index] = oa_s.AlsamSegmentStyle(
                index, row.iloc[0], traveler)

        return seg_style

    def load_parts(self):
        regions = pd.read_excel(
            self.filename, sheet_name=4, header=5, usecols="A:O", index_col=0)
        parts = pd.read_excel(self.filename, sheet_name=5,
                              header=5, usecols="A:G")

        regions_d = {}
        for index, row in regions.iterrows():
            p_r = oa_s.PrintRegions(index, *row.iloc[0:13].values.tolist())
            regions_d[index] = p_r

        part_d = []
        for index, row in parts.iterrows():
            part = oa_s.PartsSetting(*row.iloc[0:7].values.tolist())
            part_d.append(part)

        return regions_d, part_d


class PrintFileWriter():

    def __init__(self, output_folder: str, print_loader: PrintFileLoader):
        self.out = output_folder
        self.pl = print_loader

    def make_header(self):
        h = etree.Element('Header')
        etree.SubElement(h, 'AmericaMakesSchemaVersion').text = \
            self.pl.header.AmericaMakesSchemaVersion.strftime("%Y-%m-%d")
        etree.SubElement(h, 'LayerThickness').text = \
            str(self.pl.header.LayerThickness)

        return h

    def make_velocity_profile(self):
        vpl = etree.Element('VelocityProfileList')

        ls_vpl = self.pl.vel_profile

        for vl in ls_vpl.values():
            xvl = etree.SubElement(vpl, 'VelocityProfile')

            etree.SubElement(xvl, 'ID').text = str(vl.ID)
            etree.SubElement(xvl, 'Velocity').text = str(vl.Velocity)
            etree.SubElement(xvl, 'Mode').text = str(vl.Mode)
            etree.SubElement(xvl, 'LaserOnDelay').text = str(vl.LaserOnDelay)
            etree.SubElement(xvl, 'LaserOffDelay').text = str(vl.LaserOffDelay)
            etree.SubElement(xvl, 'MarkDelay').text = str(vl.MarkDelay)
            etree.SubElement(xvl, 'PolygonDelay').text = str(vl.PolygonDelay)

        return vpl

    def make_segment_style(self):
        xsl_root = etree.Element('SegmentStyleList')

        ls_seg = self.pl.seg_style
        for sl in ls_seg.values():
            xsl = etree.SubElement(xsl_root, 'SegmentStyle')

            etree.SubElement(xsl, 'ID').text = str(sl.ID)
            etree.SubElement(xsl, 'VelocityProfileID').text = str(
                sl.VelocityProfileID)
            etree.SubElement(xsl, 'LaserMode').text = str(sl.LaserMode)

            xtl = etree.SubElement(xsl_root, 'Traveler')

            etree.SubElement(xtl, 'ID').text = str(sl.Traveler.ID)
            etree.SubElement(xtl, 'SyncDelay').text = str(
                sl.Traveler.SyncDelay)
            etree.SubElement(xtl, 'Power').text = str(sl.Traveler.Power)
            etree.SubElement(xtl, 'SpotSize').text = str(sl.Traveler.SpotSize)

        return xsl_root

    def write_svg(self):
        pass

    def write(self, contour : List[Polygon], hatch : List[MultiLineString]):
        with etree.xmlfile(self.out) as xf:
            with xf.element('Layer'):
                xf.write(self.make_header())
                xf.write(self.make_velocity_profile())
                xf.write(self.make_segment_style())
