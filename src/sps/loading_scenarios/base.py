#!/usr/bin/env python3

import numpy as np
from abc import abstractmethod
from dataclasses import dataclass

@dataclass(slots=True)
class SPSOutputs:
    self.stress: list
    self.strain: list
    self.time: list
    self.Eeff: list
    self.all_dfrd: list
    self.vals: list

@dataclass(slots=True)
class BaseScenario:

