# --- STL Imports ---
import sys

if sys.argv[1:]:
    target = sys.argv[1]
else:
    target = ""

if target == "heterogeneous":
    from .heterogeneous_waveform_relaxation import run
elif target == "homogeneous":
    from .homogeneous_waveform_relaxation import run
elif target == "profile":
    from .profile import run
else:
    def run():
        pass

args = sys.argv[1:]
run(args)
