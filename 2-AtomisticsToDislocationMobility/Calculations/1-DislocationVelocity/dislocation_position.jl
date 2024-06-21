# if running WSL, then run from localhost NOT from within WSL

using Printf
using PyCall
ovito = pyimport("ovito") # must be installed on localhost Python kernel

filepath = pwd()

REFERENCE_STRUCTURE = "bcc"

# must be same examined in `rescale_commands.sh`
TEMP = [300] # range(150, 500; step=50) # K
SIGMA = range(25, 300; step=25) # MPa
# for s in [700, 800, 900, 1000, 1100, 1200]
#     push!(SIGMA, s)
# end

reference_structure = lowercase(REFERENCE_STRUCTURE)
workingdir = joinpath([filepath, "PositionFrameData"])
!(isdir(workingdir)) ? mkdir(workingdir) : nothing

for temp in TEMP
    dir_temp = joinpath([workingdir, "$temp"])
    !(isdir(dir_temp)) ? mkdir(dir_temp) : nothing
    for sigma in SIGMA
        dirtemp_sigma = joinpath([dir_temp, "$sigma"])
        !(isdir(dirtemp_sigma)) ? mkdir(dirtemp_sigma) : nothing
        @printf("Processing %3.3f/%3.3f...\r", temp, sigma)
        pipeline = ovito.io.import_file("$filepath/RescaleDownload/$temp/$sigma/dump.shear.unwrap")
        if reference_structure == "bcc"
            pipeline.modifiers.append(ovito.modifiers.CentroSymmetryModifier(num_neighbors=8))
        elseif reference_structure == "fcc"
            pipeline.modifiers.append(ovito.modifiers.CentroSymmetryModifier(num_neighbors=12))
        else
            error("Variable 'REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE' not understood. Must be either 'fcc' or 'bcc'.")
        end
        pipeline.modifiers.append(ovito.modifiers.ColorCodingModifier(property="Centrosymmetry"))
        pipeline.modifiers.append(ovito.modifiers.PolyhedralTemplateMatchingModifier())
        pipeline.modifiers.append(ovito.modifiers.ExpressionSelectionModifier(expression="(StructureType == 3 || StructureType == 2 || StructureType == 1) || (Position.Y < 20 || Position.Y > 158)"))
        pipeline.modifiers.append(ovito.modifiers.DeleteSelectedModifier())
        pipeline.modifiers.append(ovito.modifiers.ComputePropertyModifier(
            output_property="Time",
            expressions="Timestep"))
        for frame in range(0, pipeline.source.num_frames - 1)
            data = pipeline.compute(frame)
            ovito.io.export_file(data, "$filepath/PositionFrameData/$temp/$sigma/lammps.$frame.xyz", "xyz", multiple_frames=true, columns=["Particle Identifier", "Position", "Centrosymmetry", "Color", "Structure Type", "Time"])
        end
    end
end
print("\e[2K"); println("Done.")