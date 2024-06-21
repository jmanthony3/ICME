# if running WSL, then run from localhost NOT from within WSL

import numpy as np
import os
import ovito # must be installed on localhost Python kernel
import sys

filepath = os.path.dirname(os.path.abspath(__file__))

REFERENCE_STRUCTURE = "bcc"

# must be same examined in `rescale_commands.sh`
TEMP = np.array([300]) # np.arange(150, 550, 50) # K
SIGMA = np.arange(25, 325, 25) # MPa
# for s in np.array([700, 800, 900, 1000, 1100, 1200]):
#     SIGMA = np.append(SIGMA, s)

reference_structure = REFERENCE_STRUCTURE.lower()

for temp in TEMP:
    try: os.mkdir(f"{filepath}/PositionFrameData/{temp}")
    except: OSError
    for sigma in SIGMA:
        try: os.mkdir(f"{filepath}/PositionFrameData/{temp}/{sigma}")
        except: OSError
        print(f"Processing {temp}/{sigma}...", end="\r")
        pipeline = ovito.io.import_file(f"{filepath}/RescaleDownload/{temp}/{sigma}/dump.shear.unwrap")
        if reference_structure == "bcc":
            pipeline.modifiers.append(ovito.modifiers.CentroSymmetryModifier(
                num_neighbors=8
            ))
        elif reference_structure == "fcc":
            pipeline.modifiers.append(ovito.modifiers.CentroSymmetryModifier(
                num_neighbors=12
            ))
        else:
            sys.exit(f"Variable 'REFERENCE_STRUCTURE={REFERENCE_STRUCTURE}' not understood. Must be either 'fcc' or 'bcc'.")
        pipeline.modifiers.append(ovito.modifiers.ColorCodingModifier(
            property="Centrosymmetry"
        ))
        pipeline.modifiers.append(ovito.modifiers.PolyhedralTemplateMatchingModifier(
        ))
        pipeline.modifiers.append(ovito.modifiers.ExpressionSelectionModifier(
            expression="(StructureType == 3 || StructureType == 2 || StructureType == 1) || (Position.Y < 20 || Position.Y > 158)"
        ))
        pipeline.modifiers.append(ovito.modifiers.DeleteSelectedModifier(
        ))
        pipeline.modifiers.append(ovito.modifiers.ComputePropertyModifier(
            output_property="Time",
            expressions="Timestep"
        ))
        for frame in range(pipeline.source.num_frames):
            data = pipeline.compute(frame)
            ovito.io.export_file(data, f"{filepath}/PositionFrameData/{temp}/{sigma}/lammps.{frame}.xyz", "xyz", multiple_frames=True, columns=["Particle Identifier", "Position", "Centrosymmetry", "Color", "Structure Type", "Time"])
print("Done.")