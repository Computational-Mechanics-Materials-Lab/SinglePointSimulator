from typing import Self
from dataclasses import dataclass, field
import pandas as pd
import numpy.typing as npt


@dataclass(slots=True)
class SPSOutputs:
    stress: list[float] = field(default_factory=list)
    strain: list[npt.NDArray] = field(default_factory=list)
    time: list[float] = field(default_factory=list)
    Eeff: list[float] = field(default_factory=list)
    all_dfgrd: list[npt.NDArray] = field(default_factory=list)
    vals: list[float] = field(default_factory=list)

    def to_df(self: Self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "time": self.time,
                "stress": self.stress,
                "strain": self.strain,
                "Eeff": self.Eeff,
                "all_dfgrd": self.all_dfgrd,
            }
        )
