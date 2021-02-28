from dataclasses import dataclass

from datetime import date
from typing import Union, List

@dataclass
class AlsamHeader():
    AmericaMakesSchemaVersion : date
    LayerNum : int # > 0
    LayerThickness : float # > 0 optional
    AbsoluteHeight : float = -1 # > 0 optional
    DosingFactor : float = -1 # > 0 optional
    BuildDescription : str = '' # optional

@dataclass
class AlsamVelocityProfile():
    ID : str
    Velocity : float # > 0
    LaserOnDelay : Union[float, int] # microsecond
    LaserOffDelay : Union[float, int] # microsecond
    JumpDelay : Union[float, int] # microsecond
    MarkDelay : Union[float, int] 
    PolygonDelay : Union[float, int]
    Mode : str = 'Auto' # Delay or Auto

@dataclass
class AlsamWobble():
    On : int # 0 or 1
    Freq: int 
    Shape: int # [-1, 0, 1]
    TransAmp : float 
    LongAmp : float  

@dataclass
class AlsamTraveler():
    ID : str 
    SyncDelay : int # only if multi laser
    Power : float = -1 # optional
    SpotSize : float = -1 # optional
    Wobble : AlsamWobble = None # optional

@dataclass
class AlsamSegmentStyle():
    ID : str
    VelocityProfileID : str 
    LaserMode : str = "Independent" # FollowMe or Independent
    Traveler : AlsamTraveler

@dataclass
class AlsamSegStart():
    X : float
    Y : float

@dataclass
class AlsamSegEnd():
    X : float 
    Y : float

@dataclass
class AlsamSegment():
    SegmentID : str
    SegStyle : str
    End : AlsamSegEnd

@dataclass
class AlsamPath():
    Type : str # hatch or contour
    Tag : str 
    NumSegments : int 
    SkyWritingMode : int # 0, 1, 2, 3
    Start : AlsamSegStart
    Segment : List[AlsamSegment]

@dataclass 
class AlsamTrajectory():
    TrajectoryID : int
    PathProcessingMode : str # sequential or concurrent
    Path : List[AlsamPath]