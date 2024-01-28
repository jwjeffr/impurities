#!/usr/bin/python

"""
Script for creating MC-MD figure
"""

from functools import partial
import ovito
from PIL import Image
from constants import AtomColors, Alignments, SYSTEMS


def type_map_modifier(
    frame: int, data: ovito.data.DataCollection, type_info_map: dict
) -> None:
    """
    modifier that adds particle names and radii to a dump file pipeline
    important for labeling the type legend
    :param frame: frame of pipeline
    :param data: data collection from frame
    :param type_info_map: dictionary mapping integer labels to (atom type name, radius, display color) tuples
    :return: None
    """

    # get integer labels of types
    types = data.particles_.particle_types_

    # assign atom type name and radius to types
    for key, (atom_type_name, radius, color) in type_info_map.items():
        types.type_by_id_(key).name = atom_type_name
        types.type_by_id_(key).radius = radius
        types.type_by_id_(key).color = color


def generate_image(system: str) -> Image:
    # arguments for saving ovito rendered images
    saving_keyword_arguments = dict(
        size=(1000, 1000), renderer=ovito.vis.TachyonRenderer(), alpha=True
    )

    # type radius map and its associated modifier
    if system == "cantor":
        type_info_map = {
            1: ("Co", 1.25, AtomColors.COBALT),
            2: ("Ni", 1.25, AtomColors.NICKEL),
            3: ("Cr", 1.29, AtomColors.CHROMIUM),
            4: ("Fe", 1.26, AtomColors.IRON),
            5: ("Mn", 1.39, AtomColors.MANGANESE),
        }
    elif system == "FeAl":
        type_info_map = {
            1: ("Fe", 1.26, AtomColors.IRON),
            2: ("Al", 1.43, AtomColors.ALUMINUM),
        }
    modifier = partial(type_map_modifier, type_info_map=type_info_map)

    for state in ["initial", "final"]:
        # initialize pipeline, add in modifier, add it to the scene for rendering
        pipeline = ovito.io.import_file(f"mc_data/{system}/mc.dump")
        pipeline.modifiers.append(modifier)
        pipeline.add_to_scene()

        # create a viewport for rendering
        vp = ovito.vis.Viewport(
            type=ovito.vis.Viewport.Type.Perspective, camera_dir=(-1, -1, -0.75)
        )

        # create either a type legend or a coordinate tripod depending on the state being rendered
        # so type legend is at the top of the image and coordinate tripod is on the bottom of the image
        if state == "initial":
            overlay = ovito.vis.ColorLegendOverlay(
                title=" ",
                alignment=Alignments.TOP_LEFT,
                offset_y=-0.06,
                offset_x=0.02,
                font_size=0.1,
                property="particles/Particle Type",
            )
        elif state == "final" and system == "cantor":
            overlay = ovito.vis.CoordinateTripodOverlay(
                alignment=Alignments.BOTTOM_LEFT,
                axis1_label="100",
                axis1_color=(0, 0, 0),
                axis2_label="010",
                axis2_color=(0, 0, 0),
                axis3_label="001",
                axis3_color=(0, 0, 0),
                offset_x=0.1,
                size=0.08,
            )
        else:
            overlay = None

        # add overlay, render image at the appropriate frame
        try:
            vp.overlays.append(overlay)
        except ValueError:
            pass
        vp.zoom_all()
        if state == "initial":
            frame = 0
        if state == "final":
            frame = pipeline.source.num_frames
        vp.render_image(
            filename=f"plots/{system}_mc_{state}.png",
            frame=frame,
            **saving_keyword_arguments,
        )

        # remove current pipeline from scene so we can add the final state
        pipeline.remove_from_scene()

    # open images to combine them
    initial_image = Image.open(f"plots/{system}_mc_initial.png")
    final_image = Image.open(f"plots/{system}_mc_final.png")

    # get widths and heights
    widths, heights = zip(*(i.size for i in [initial_image, final_image]))

    # create a crop variable, a lot of empty space above both images
    crop = 125

    # get height and width of combined image
    total_height = sum(heights) - crop
    max_width = max(widths)

    # initialize a new image, paste final state, paste initial state over the empty space above the final state's image
    new_im = Image.new("RGBA", (max_width, total_height))
    new_im.paste(final_image, (0, initial_image.size[1] - crop))
    new_im.paste(initial_image, (0, 0))

    return new_im


def main():
    images = [generate_image(system) for system in SYSTEMS]

    widths, heights = zip(*(i.size for i in images))

    # get height and width of combined image
    total_width = sum(widths)
    max_height = max(heights)

    # initialize a new image, paste final state, paste initial state over the empty space above the final state's image
    new_im = Image.new("RGBA", (total_width, max_height))
    new_im.paste(images[0], (0, 0))
    new_im.paste(images[1], (images[1].size[0], 0))

    new_im.save("plots/mc.png")


if __name__ == "__main__":
    main()
