
import sys, subprocess
from pathlib import Path

#*----------------------------------------------------------------------------*#

    ########################################################################

'''  
      revision log:

        18 Apr. 2026 (Hao Li)
          --- initial commit: implement fits batch inversion
'''

    ########################################################################

#*----------------------------------------------------------------------------*#

if __name__ == "__main__":
    
    if len(sys.argv) < 3:
      raise ValueError("Usage: python RunMEINV.py <data_path> np=<core_number>")

    path = sys.argv[1].strip()
    path = Path(path)

    print('path = ', path)
    if not path.exists():
      raise ValueError("path does not exist")

    nproc = 6
    if 'np=' in sys.argv[2].lower():
      try:
        nproc = int(sys.argv[2].split("=")[1])
      except:
        raise ValueError("Invalid np format")
 
    print("nproc =", nproc)     


    files = sorted(path.glob("*.fits"))

    if len(files) == 0:
      raise ValueError("No data file found")

    outdir = path / "res"
    outdir.mkdir(exist_ok=True)

    for i, f in enumerate(files):

      stem = f.stem
      output = outdir / f"{stem}_res.fits"

      config = {
        "verbose": 0,
        "num_run": 5,
        "lines": "15648.5, 3.0",
        "wavelength_path": "../niris_example/niriswavelength.fits",
        "weights": "1, 8, 8, 12",
        "Bounds_Bmod": "5, 5000",
        "Bounds_Vlos": "-20, 20",
        "Bounds_Dopp": "50, 140",
        "Bounds_Damp": "0.4, 1.5",
        "Bounds_Eta": "2, 50",
        "Bounds_Beta": "0.1, 0.9",
        "chisq_criteria": 5,
        "LCoeffi": 330,
        "VCoeffi": 330,
        "result_path": str(output),
        "data_path": str(f)
      }

      with open("input", "w") as fp:
        for k, v in config.items():
          fp.write(f"{k} = {v}\n")

      print(f"Processing: {i+1}/{len(files)}  {f.name}")

      cmd = ["mpirun", "-np", str(nproc), "../MEINV"]
      subprocess.run(cmd, check=True)
  
