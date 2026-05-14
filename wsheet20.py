"""
wsheet20.py  —  Python wrapper for WSHEET20.F
Wide-band Thin Sheet EM Integral Equation modelling (whole space).

Usage
-----
    from wsheet20 import Wsheet20, Sheet, ColeCole, Source, Receivers

    run = Wsheet20(
        sheets=[
            Sheet(
                na=20, nb=20,
                center=(500., 0., -200.),
                strike_length=200., dip_length=200.,
                strike=90., dip=90., thickness=0.1,
                resistivity=ColeCole(low=1e-5, high=1e-6, tau=0., alpha=0.5),
                permittivity=ColeCole(low=2., high=2., tau=0., alpha=0.5),
            )
        ],
        host_resistivity=ColeCole(low=10., high=10., tau=0., alpha=0.5),
        host_permittivity=ColeCole(low=10., high=10., tau=0., alpha=0.5),
        frequency=20e3,
        source=Source(type="magnetic_dipole", direction="z", angle=0.),
        tx_start=(0., 0., 0.), tx_increment=(0., 0., 0.), n_tx=1,
        rx_start=(0., 0., 0.75), rx_increment=(5., 0., 0.), n_rx=201,
        green_mode=2,           # 0=compute+save, 1=compute only, 2=reuse saved
        output_scatter=True,    # write scatter/incident fields
        output_efield=True,     # write e-field.dat
        output_hfield=False,    # write h-field.dat
        normalize=False,        # normalise by primary field
    )
    run.write_par()             # writes wsheet20.par in the working directory
    result = run.run()          # runs ./a.out and returns CompletedProcess
"""

from __future__ import annotations

import subprocess
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Tuple


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ColeCole:
    """Cole-Cole dispersion parameters.

    Parameters
    ----------
    low :   low-frequency limit  (z_l in Fortran)
    high:   high-frequency limit (z_h)   — must be <= low
    tau :   time constant (s);  set to 0 for non-polarisable
    alpha:  exponent (0 < alpha <= 1)
    """
    low:   float
    high:  float
    tau:   float = 0.0
    alpha: float = 0.5

    def _fmt(self) -> str:
        """Fortran-style: low high tau alpha   (e.g. 1.d-5  1.d-6  0.  0.5)"""
        def _f(v: float) -> str:
            # Use 'd' exponent notation when the exponent is present
            s = f"{v:.6g}"
            if "e" in s or "E" in s:
                s = s.replace("e", "d").replace("E", "d")
            return s
        return f"{_f(self.low)}  {_f(self.high)}  {_f(self.tau)}  {_f(self.alpha)}"


@dataclass
class Sheet:
    """One thin sheet.

    Parameters
    ----------
    na, nb          : cells along strike / dip direction  (max 20 each)
    center          : (x, y, z) of sheet centre before rotation (m)
    strike_length   : along-strike dimension (m)
    dip_length      : down-dip dimension (m)
    strike          : azimuth from x-axis (degrees, 0–360)
    dip             : dip angle from horizontal (degrees, 0–90)
    thickness       : sheet thickness (m)
    resistivity     : ColeCole for the sheet
    permittivity    : ColeCole (relative) for the sheet
    """
    na:           int   = 20
    nb:           int   = 20
    center:       Tuple[float, float, float] = (0., 0., -200.)
    strike_length: float = 200.
    dip_length:   float  = 200.
    strike:       float  = 90.
    dip:          float  = 90.
    thickness:    float  = 0.1
    resistivity:  ColeCole = field(default_factory=lambda: ColeCole(1e-5, 1e-6))
    permittivity: ColeCole = field(default_factory=lambda: ColeCole(2., 2.))


@dataclass
class Source:
    """EM source parameters.

    Parameters
    ----------
    type      : "plane_wave" | "electric_dipole" | "magnetic_dipole"
    direction : for dipoles "x"/"y"/"z";
                for plane wave "TM" (maps to 1) or "TE" (maps to 2)
    angle     : incidence angle for plane wave (degrees); ignored for dipoles
    """
    type:      Literal["plane_wave", "electric_dipole", "magnetic_dipole"] = "magnetic_dipole"
    direction: Literal["x", "y", "z", "TM", "TE"] = "z"
    angle:     float = 0.

    _TYPE_MAP = {"plane_wave": 0, "electric_dipole": 1, "magnetic_dipole": 2}
    _DIR_MAP  = {"x": 1, "y": 2, "z": 3, "TM": 1, "TE": 2}

    def _type_int(self) -> int:
        return self._TYPE_MAP[self.type]

    def _dir_int(self) -> int:
        return self._DIR_MAP[self.direction]


# ---------------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------------

class Wsheet20:
    """Python interface to the WSHEET20 Fortran code.

    Parameters
    ----------
    sheets           : list of 1 or 2 Sheet objects
    host_resistivity : ColeCole for the whole-space resistivity
    host_permittivity: ColeCole for the whole-space relative permittivity
    frequency        : frequency in Hz
    source           : Source object
    n_tx             : number of transmitter positions
    tx_start         : (xs, ys, zs) — starting transmitter position (m)
    tx_increment     : (dxs, dys, dzs) — Tx position step (m)
    n_rx             : number of receiver positions
    rx_start         : (xr, yr, zr) — starting receiver position (m)
    rx_increment     : (dxr, dyr, dzr) — Rx position step (m)
    green_mode       : 0 = compute Green + write LUD
                       1 = compute Green, do NOT write LUD
                       2 = reuse previously written LUD (fastest for repeat runs)
    output_scatter   : write incident.out / scatter.out
    output_efield    : write e-field.dat
    output_hfield    : write h-field.dat
    normalize        : normalise secondary/total field by primary field
    executable       : path to the compiled binary (default: "./a.out")
    par_file         : path for the parameter file (default: "wsheet20.par")
    work_dir         : working directory for the run (default: current dir)
    """

    def __init__(
        self,
        sheets:            list[Sheet]       = None,
        host_resistivity:  ColeCole          = None,
        host_permittivity: ColeCole          = None,
        frequency:         float             = 20e3,
        source:            Source            = None,
        n_tx:              int               = 1,
        tx_start:          tuple             = (0., 0., 0.),
        tx_increment:      tuple             = (0., 0., 0.),
        n_rx:              int               = 201,
        rx_start:          tuple             = (0., 0., 0.75),
        rx_increment:      tuple             = (5., 0., 0.),
        green_mode:        int               = 2,
        output_scatter:    bool              = True,
        output_efield:     bool              = True,
        output_hfield:     bool              = False,
        normalize:         bool              = False,
        executable:        str | Path        = "./a.out",
        par_file:          str | Path        = "wsheet20.par",
        work_dir:          str | Path        = ".",
    ):
        if sheets is None:
            sheets = [Sheet()]
        if len(sheets) not in (1, 2):
            raise ValueError("WSHEET20 supports 1 or 2 sheets only.")

        self.sheets            = sheets
        self.host_resistivity  = host_resistivity  or ColeCole(10., 10.)
        self.host_permittivity = host_permittivity or ColeCole(10., 10.)
        self.frequency         = frequency
        self.source            = source or Source()
        self.n_tx              = n_tx
        self.tx_start          = tuple(tx_start)
        self.tx_increment      = tuple(tx_increment)
        self.n_rx              = n_rx
        self.rx_start          = tuple(rx_start)
        self.rx_increment      = tuple(rx_increment)
        self.green_mode        = green_mode
        self.output_scatter    = output_scatter
        self.output_efield     = output_efield
        self.output_hfield     = output_hfield
        self.normalize         = normalize
        self.executable        = Path(executable)
        self.par_file          = Path(par_file)
        self.work_dir          = Path(work_dir)

    # ------------------------------------------------------------------
    # Parameter-file generation
    # ------------------------------------------------------------------

    def _build_par(self) -> str:
        ns = len(self.sheets)
        lines = []

        def ln(comment, *values):
            lines.append(f">> {comment}")
            lines.append("  ".join(str(v) for v in values))

        # No. of sheets
        ln("No. of sheet", ns)

        # Cell counts
        lines.append(">> No. of cells for strike and dip direction (max = 20)")
        for s in self.sheets:
            lines.append(f"{s.na}  {s.nb}")

        # Sheet centres
        lines.append(">> xyz coordinate of the center of the sheets before rotation")
        for s in self.sheets:
            x, y, z = s.center
            lines.append(f"{x:.6g}  {y:.6g}  {z:.6g}")

        # Geometry
        lines.append(">> strike length (m), dip length(m), strike (deg), dip (deg),thickness(m) of each sheet")
        for s in self.sheets:
            lines.append(
                f"{s.strike_length:.6g}  {s.dip_length:.6g}  "
                f"{s.strike:.6g}  {s.dip:.6g}  {s.thickness:.6g}"
            )

        # Cole-Cole resistivity of sheets
        lines.append(">> Cole-Cole parameter of resistivity of each sheet")
        for s in self.sheets:
            lines.append(s.resistivity._fmt())

        # Cole-Cole permittivity of sheets
        lines.append(">> Cole-Cole parameter of relative permittivity of each sheet")
        for s in self.sheets:
            lines.append(s.permittivity._fmt())

        # Host
        lines.append(">> Cole-Cole parameter of resistivity of whole space")
        lines.append(self.host_resistivity._fmt())
        lines.append(">> Cole-Cole parameter of dielectric constant of whole space")
        lines.append(self.host_permittivity._fmt())

        # Frequency
        lines.append(">> Frequency (Hz)")
        lines.append(f"{self.frequency:.6g}")

        # Source
        src = self.source
        lines.append(">> source type (0=plane, 1=J type, 2=M type), direction (1=x,2=y,3=z), angle")
        lines.append(f"{src._type_int()}   {src._dir_int()}   {src.angle:.6g}")

        # Tx
        ln("No. of Tx", self.n_tx)
        lines.append(">> start x,y,z coordinates of the sources and increments")
        lines.append("  ".join(f"{v:.6g}" for v in self.tx_start))
        lines.append("  ".join(f"{v:.6g}" for v in self.tx_increment))

        # Rx
        ln("No. of Rx", self.n_rx)
        lines.append(">> start x,y,z coordinates of the Receivers and increments")
        lines.append("  ".join(f"{v:.6g}" for v in self.rx_start))
        lines.append("  ".join(f"{v:.6g}" for v in self.rx_increment))

        # Green / LUD control
        lines.append(">> Perform Green's function integral and SVD (0,1) or use previous SVD data(2)")
        lines.append(str(self.green_mode))

        # Output control  (nprt 1-4)
        # nprt(1): 0=write scatter/incident, 1=skip
        # nprt(2): 0=write e-field,          1=skip
        # nprt(3): 0=write h-field,          1=skip
        # nprt(4): 0=absolute,               1=normalised
        nprt1 = 0 if self.output_scatter else 1
        nprt2 = 0 if self.output_efield  else 1
        nprt3 = 0 if self.output_hfield  else 1
        nprt4 = 1 if self.normalize      else 0
        lines.append(">> output control parameter scattering, e_field, h_field, normalizing")
        lines.append(f"{nprt1}  {nprt2}  {nprt3}  {nprt4}")

        return "\n".join(lines) + "\n"

    def write_par(self, path: str | Path = None) -> Path:
        """Write the parameter file and return the path used."""
        dest = Path(path) if path else self.work_dir / self.par_file
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(self._build_par())
        print(f"[wsheet20] Parameter file written → {dest}")
        return dest

    # ------------------------------------------------------------------
    # Running the binary
    # ------------------------------------------------------------------

    def run(
        self,
        extra_args: list[str] = None,
        capture_output: bool  = False,
        timeout: float        = None,
    ) -> subprocess.CompletedProcess:
        """Write the .par file and execute the Fortran binary.

        Parameters
        ----------
        extra_args     : additional CLI arguments passed to the binary
        capture_output : if True, stdout/stderr are captured (not printed)
        timeout        : optional timeout in seconds

        Returns
        -------
        subprocess.CompletedProcess
        """
        par_path = self.write_par()

        exe = self.executable
        if not exe.is_absolute():
            exe = (self.work_dir / exe).resolve()

        if not exe.exists():
            raise FileNotFoundError(
                f"Executable not found: {exe}\n"
                "Compile with e.g.:  gfortran wsheet20.f SUB_LUD.f -o a.out"
            )

        cmd = [str(exe)] + (extra_args or [])
        print(f"[wsheet20] Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            cwd=str(self.work_dir),
            capture_output=capture_output,
            text=True,
            timeout=timeout,
        )

        if result.returncode != 0:
            print(f"[wsheet20] *** Non-zero exit code: {result.returncode} ***")
        else:
            print("[wsheet20] Finished successfully.")

        return result

    # ------------------------------------------------------------------
    # Convenience: show current parameters
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a human-readable summary of all parameters."""
        src = self.source
        lines = [
            "=" * 60,
            "WSHEET20 run parameters",
            "=" * 60,
            f"  Sheets          : {len(self.sheets)}",
        ]
        for i, s in enumerate(self.sheets, 1):
            lines += [
                f"  --- Sheet {i} ---",
                f"    Grid          : {s.na} × {s.nb} cells",
                f"    Centre (m)    : {s.center}",
                f"    Strike len.   : {s.strike_length} m",
                f"    Dip len.      : {s.dip_length} m",
                f"    Strike / Dip  : {s.strike}° / {s.dip}°",
                f"    Thickness     : {s.thickness} m",
                f"    Resistivity   : {s.resistivity}",
                f"    Permittivity  : {s.permittivity}",
            ]
        lines += [
            f"  Host resist.    : {self.host_resistivity}",
            f"  Host permitt.   : {self.host_permittivity}",
            f"  Frequency       : {self.frequency:.6g} Hz",
            f"  Source type     : {src.type}",
            f"  Source dir.     : {src.direction}  (angle={src.angle}°)",
            f"  n_tx / n_rx     : {self.n_tx} / {self.n_rx}",
            f"  Tx start        : {self.tx_start}  Δ={self.tx_increment}",
            f"  Rx start        : {self.rx_start}  Δ={self.rx_increment}",
            f"  Green mode      : {self.green_mode}",
            f"  Output scatter  : {self.output_scatter}",
            f"  Output E-field  : {self.output_efield}",
            f"  Output H-field  : {self.output_hfield}",
            f"  Normalise       : {self.normalize}",
            "=" * 60,
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()


# ---------------------------------------------------------------------------
# Quick-start example (matches the original wsheet20.par)
# ---------------------------------------------------------------------------

def _default_run() -> Wsheet20:
    """Recreate the original wsheet20.par settings."""
    return Wsheet20(
        sheets=[
            Sheet(
                na=20, nb=20,
                center=(500., 0., -200.),
                strike_length=200., dip_length=200.,
                strike=90., dip=90., thickness=0.1,
                resistivity=ColeCole(low=1e-5, high=1e-6, tau=0., alpha=0.5),
                permittivity=ColeCole(low=2., high=2., tau=0., alpha=0.5),
            )
        ],
        host_resistivity  = ColeCole(low=10., high=10., tau=0., alpha=0.5),
        host_permittivity = ColeCole(low=10., high=10., tau=0., alpha=0.5),
        frequency    = 20e3,
        source       = Source(type="magnetic_dipole", direction="z", angle=0.),
        n_tx         = 1,
        tx_start     = (0., 0., 0.),
        tx_increment = (0., 0., 0.),
        n_rx         = 201,
        rx_start     = (0., 0., 0.75),
        rx_increment = (5., 0., 0.),
        green_mode   = 2,
        output_scatter = True,
        output_efield  = True,
        output_hfield  = False,
        normalize      = False,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Python wrapper for WSHEET20.F EM thin-sheet modelling code."
    )
    parser.add_argument("--write-par-only", action="store_true",
                        help="Write parameter file and exit without running.")
    parser.add_argument("--exe", default="./a.out",
                        help="Path to the compiled Fortran binary (default: ./a.out)")
    parser.add_argument("--work-dir", default=".",
                        help="Working directory (default: current)")
    args = parser.parse_args()

    run = _default_run()
    run.executable = Path(args.exe)
    run.work_dir   = Path(args.work_dir)

    print(run.summary())

    if args.write_par_only:
        run.write_par()
    else:
        run.run()
