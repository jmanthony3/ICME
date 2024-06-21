# if running WSL, then run from localhost NOT from within WSL

using Printf
using PyCall
ovito = pyimport("ovito") # must be installed on localhost Python kernel

filepath = pwd()

REFERENCE_STRUCTURE = "bcc"

# must be same examined in `rescale_commands.sh`
TEMP = [300.] # float.(range(150, 500; step=50)) # K
SIGMA = float.(range(25, 300; step=25)) # MPa
# append!(SIGMA, float.([700, 800, 900, 1000, 1100, 1200]))

reference_structure = lowercase(REFERENCE_STRUCTURE)
dir_pfd = joinpath([filepath, "PositionFrameData"])
# dir_rescale_up = joinpath([filepath, "RescaleUpload"])
dir_rescale_down = joinpath([filepath, "RescaleDownload"])
!(isdir(dir_pfd)) ? mkdir(dir_pfd) : nothing
# !(isdir(dir_rescale_up)) ? mkdir(dir_rescale_up) : nothing
!(isdir(dir_rescale_down)) ? mkdir(dir_rescale_down) : nothing

for temp in TEMP
    dirpfd_temp = joinpath([dir_pfd, @sprintf("%.3f", temp)])
    # dirrscl_up = joinpath([dir_rescale_up, @sprintf("%.3f", temp)])
    dirrscl_down = joinpath([dir_rescale_down, @sprintf("%.3f", temp)])
    !(isdir(dirpfd_temp)) ? mkdir(dirpfd_temp) : nothing
    # !(isdir(dirrscl_up)) ? mkdir(dirrscl_up) : nothing
    !(isdir(dirrscl_down)) ? mkdir(dirrscl_down) : nothing
    for sigma in SIGMA
        dirpfdtemp_sigma = joinpath([dirpfd_temp, @sprintf("%.3f", sigma)])
        # diruptemp_sigma = joinpath([dirrscl_up, @sprintf("%.3f", sigma)])
        dirdowntemp_sigma = joinpath([dirrscl_down, @sprintf("%.3f", sigma)])
        !(isdir(dirpfdtemp_sigma)) ? mkdir(dirpfdtemp_sigma) : nothing
        # !(isdir(diruptemp_sigma)) ? mkdir(diruptemp_sigma) : nothing
        !(isdir(dirdowntemp_sigma)) ? mkdir(dirdowntemp_sigma) : nothing
        @printf("Processing %.3f/%.3f...\r", temp, sigma)
        # println(dirdowntemp_sigma)
        pipeline = ovito.io.import_file("$dirdowntemp_sigma/dump.shear.unwrap")
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
            ovito.io.export_file(data, "$dirpfdtemp_sigma/lammps.$frame.xyz", "xyz", multiple_frames=true, columns=["Particle Identifier", "Position", "Centrosymmetry", "Color", "Structure Type", "Time"])
        end
    end
end
print("\e[2K"); println("Done.")