import os
import Sofa
import SofaRuntime
from SofaRuntime import Timer
import time

timers = {}

def main(root_dir="."):
    for root, dirs, files in os.walk(root_dir):
        if "beam.scn" in files:
            beam_path = os.path.join(root, "beam.scn")

            # Extra safety check
            if not os.path.isfile(beam_path):
                continue

            print(f"Running {beam_path}")

            root = Sofa.Simulation.load(beam_path)

            # Once defined, initialization of the scene graph
            Sofa.Simulation.initRoot(root)

            Timer.setEnabled("Animate", True)
            Timer.setOutputType("Animate", "json")  # hack to avoid printing in the console

            start = time.perf_counter()

            # Run as many simulation steps
            for iteration in range(100):
                Timer.begin("Animate")
                Sofa.Simulation.animate(root, root.dt.value)
                Timer.end("Animate")

            end = time.perf_counter()
            elapsed = end - start

            print(f"Elapsed {elapsed} s")

            timers[beam_path] = elapsed

            Timer.clear()
            Sofa.Simulation.unload(root)

if __name__ == "__main__":
    main()

    for scene,t in timers.items():
        print(f"{scene}: {t}")
