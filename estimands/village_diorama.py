"""Stylised village diorama with treatment and spillover for Blender 5.1.

In Blender: Scripting > Text > Open this file > Run Script (Alt+P).
Then Render > Render Image (F12). Save the project with File > Save As.

Everything is generated locally: geometry, materials, lighting and camera.
No add-ons, image textures, internet access or external Python packages.
Each run creates a NEW scene; existing scenes and objects are preserved.
The model is schematic: house sizes and village circles are illustrative.
Distances shown are village-centre distances in arbitrary scene units.
Blue roofs and flags mean direct treatment. Blue shells, terraces, contours,
mist or ground tint mean illustrative spillover intensity S(x,y), not direct
assignment or a fitted outcome.

Spillover display styles (SPILLOVER_STYLE, combine any of them with "+"):
  "shells"   Nested translucent domes, one per value in FIELD_LEVELS. Each dome
             is an iso-surface of S(x,y) * exp(-(z/SHELL_HEIGHT)^2), so its
             footprint is exactly the ground contour S(x,y) = level. Discrete
             nested steps read far more clearly than a continuous haze.
  "terraces" Stepped glass: a wall on each contour, taller for inner levels,
             over a floor tint that steps up towards the sources. Less
             volumetric than shells, but never veils the houses from above.
  "contours" A coloured ground line along each contour, so every house is
             unambiguously inside or outside each level.
  "mist"     The original continuous volume. Its opacity integrates density
             along the viewing ray, which smears the gradient sideways, so it
             is not a pointwise effect scale.
  "ground"   A continuous planar colour gradient on the ground.
Shells, terraces and contours use exact contours and are cut off at the board
edge where a level extends beyond it, as if the landscape continued.
Keep geography, kernel, source weights, levels and colours fixed across panels.
Change SCENARIO to switch counterfactual assignments on the same landscape.

Lighting and output: a warm Sun lamp key, cool fill and cool rim light, soft
skylight from Blender's built-in Sky Texture, and an AgX contrast look. The
camera zoom is computed from the board alone, so every scenario and style gets
the identical frame and panels can be stacked or animated without shifting.
BACKGROUND = "transparent" writes a PNG with alpha (the board's shadow kept) for
placing on a poster or page colour. Labels are off by default: the console lists
pixel positions of villages, houses and sources for adding text in the document.
"""

import json
import math
import random
import bpy
import bmesh
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector, Matrix


# -------------------------- Editable settings --------------------------
SEED = 29
RESOLUTION = (1500, 1150)
SAMPLES = 64                       # Try 24 for a quicker draft, 128 for final.
RENDER_NOW = False                 # Otherwise press F12 when ready.
SHOW_VILLAGE_RINGS = True
SHOW_TREATMENT_FLAGS = True        # A shape cue in addition to blue roof colour.
SHOW_SPILLOVER = True
SPILLOVER_STYLE = "shells+contours"  # Any "+" combination of "shells",
                                     # "terraces", "contours", "mist", "ground",
                                     # e.g. "terraces+contours".
                                     # "both" is kept as "mist+ground".
SHOW_LABELS = False                # Text is best added in the document itself; the
                                   # console prints pixel anchors for placing it.
SHOW_DISTANCE_GUIDES = True       # Dashed village-centre lines (numbers need labels).
VILLAGE_RADIUS = 2.15
GROUND_SIZE = (14.8, 12.0)

# Lighting, colour and framing. All fixed across panels.
BACKGROUND = "studio"              # "studio": cool grey backdrop with a vignette.
                                   # "transparent": PNG with alpha that keeps the
                                   # board's shadow, to drop onto any page colour.
BACKDROP_COLOURS = ("C3C8CE", "8A929A")  # Backdrop centre -> frame corners.
SUN_STRENGTH = 3.2                 # Warm key: a Sun lamp, for crisp shadows.
SUN_ANGLE = 3.0                    # Degrees; smaller = sharper shadow edges.
SUN_POSITION = (-7, -9, 14)        # Direction the sun shines from (front-left, high).
FILL_RATIO = .33                   # Cool fill from the camera's right, vs the key.
RIM_STRENGTH = .8                 # Cool light from behind: catches roof ridges and
                                   # tree crowns, lifts the board off the backdrop.
SKY_STRENGTH = .07                 # Built-in Sky Texture as soft skylight only.
COLOUR_LOOK = "AgX - Medium High Contrast"  # Or "AgX - Base Contrast" (previous),
                                   # "AgX - High Contrast", "AgX - Punchy".
EXPOSURE = .1
FRAME_MARGIN = .035                # Space left around the board, per side, as a
                                   # fraction of the frame. The camera zoom and
                                   # centring are computed from the board, never from
                                   # treatment or spillover, so panels line up exactly.
LABEL_ANCHORS_FILE = ""            # e.g. "//anchors.json" (next to the saved .blend):
                                   # pixel positions of villages, houses and sources.

# Landscape positions; z is height, so these are x,y coordinates on the ground.
# A is in the foreground. A--B = 4.776; A--C = 9.192 scene units.
VILLAGES = [
    {"name": "A", "xy": ( 2.1, -3.4), "role": "focal"},
    {"name": "B", "xy": ( 3.7,  1.1), "role": "near"},
    {"name": "C", "xy": (-4.4,  3.1), "role": "far"},
]

# Each entry gives x, y, rotation in degrees and scale within a village.
# House numbering starts at 1. Add/remove entries to change the house count.
HOUSE_LAYOUT = [(-.97,.40,-9,.88), (.90,.53,12,.83), (-.10,-.99,-5,.81)]

SCENARIO = "near_village"
# Options: "untreated", "single_house", "two_houses", "near_village",
#          "far_village", "separate_source", "manual".

# For SCENARIO = "separate_source": a non-house source (e.g. a water point or
# clinic) standing INSIDE one village circle. Choose "A", "B" or "C".
# The offset is relative to that village centre, in the open ground in front of
# the houses, facing the camera. The same offset is used in every village, so the
# only thing that changes between panels is which village hosts the source.
SEPARATE_SOURCE_VILLAGE = "B"
SEPARATE_SOURCE_OFFSET = (1.00, -1.30)
SEPARATE_SOURCE_STRENGTH = 1.0

# For SCENARIO = "manual", use any selection, e.g. {"A": [1], "B": [1, 3]}.
# Use "all" to treat every house in the named village.
MANUAL_TREATED_HOUSES = {"B": [1, 3]}
MANUAL_EXTERNAL_SOURCES = [
    # {"name": "S1", "village": "B", "offset": (-2.7, 0), "strength": 1.0},
    # Or an absolute position: {"name": "S1", "xy": (1.0, 1.1), "strength": 1.0}
]

# "houses": each treated house is one additive source of HOUSE_SOURCE_STRENGTH.
# "village_centres": each fully treated village is ONE centre source of
# VILLAGE_SOURCE_STRENGTH; partial treatment is then rejected as ambiguous.
# Thus treating three houses is not silently equated with one unit source.
SOURCE_UNIT = "houses"
HOUSE_SOURCE_STRENGTH = 1.0
VILLAGE_SOURCE_STRENGTH = 1.0
KERNEL = "exponential"             # Or "gaussian".
DECAY_LENGTH = 2.2                  # Exponential: intensity falls to 1/e here.
FIELD_COLOUR = "389BE6"             # Change to e.g. "DE6653" for warm red.
MAX_FIELD_OPACITY = .78
FIELD_REFERENCE = .8
# Field S(x) = sum strength_j * kernel(distance(x, source_j), length_j).
# Display opacity = MAX_FIELD_OPACITY * (1-exp(-S(x)/FIELD_REFERENCE)).
# This is a fixed, monotone display mapping, NOT panel-specific normalisation.
# The smooth mapping compresses high intensities; the picture is illustrative.

# Shell and contour controls. FIELD_LEVELS are ABSOLUTE values of S(x,y), fixed
# across panels, innermost first. Halving steps make each shell read as "half
# the intensity of the one inside it". A single unit-strength house source never
# exceeds 1, so its innermost shell is small; three treated houses exceed .8
# across most of their village.
FIELD_LEVELS = (.8, .4, .2, .1)
FIELD_LEVEL_COLOURS = ("1F5FAE", "2F80CF", "5BA6E3", "93C8F0")  # Inner -> outer.
SHELL_HEIGHT = .95                 # Dome height scale: z = H*sqrt(ln(S/level)).
SHELL_OPACITY = .07                # Opacity where a dome faces the camera...
SHELL_RIM = .55                    # ...plus this much at grazing angles, so the
                                   # silhouette of each dome reads as a surface.
SHELL_GLOW = .5                   # Emission keeps the blue from greying out.
TERRACE_STEP = .30                 # "terraces": wall height added per level inward.
TERRACE_OPACITY = .45              # Wall opacity at its foot, thinning upward.
TERRACE_FLOOR = .10                # Tint added on the ground by each level, so
                                   # overlapping floors step up towards sources.
SHELL_GRID = .05                   # Mesh resolution of domes/contours (scene units).
CONTOUR_RADIUS = .024              # Thickness of the ground contour lines.

# Mist controls: spatial decay is still governed by DECAY_LENGTH above.
# Increase density to .40 for stronger fog; reduce to .15 for more transparent fog.
MIST_DENSITY = .28                 # Global multiplier of the SUM of source kernels.
MIST_HEIGHT = 1.4                  # Vertical e-folding height, in scene units.
MIST_COLOUR = "69B9EE"             # Luminous blue; independent of treatment roofs.
MIST_GLOW = .30                    # Small emission component makes the blue visible.
MIST_EDGE_FADE = .65               # Smooth DISPLAY crop at the landscape edges.
# density(x,y,z) = MIST_DENSITY * S(x,y) * exp(-(z/MIST_HEIGHT)^4)
#                 * edge_crop(x,y).
# The same vertical envelope, crop and density scale are used in every panel.
# No random turbulence: additional spatial variation would obscure the kernel.

SPILLOVER_STYLES = set(("mist+ground" if SPILLOVER_STYLE == "both" else SPILLOVER_STYLE)
                       .replace(" ", "").split("+")) - {""}

# Fixed village furniture, in village-local coordinates (shared with the builder
# below so the source-placement check always matches what is drawn).
VILLAGE_TREE = ((-.30, 1.35), .90)                     # position, scale
VILLAGE_SHRUBS = [((1.80*math.cos(.12+j*.70), 1.80*math.sin(.12+j*.70)), .24)
                  for j in range(4)]                   # position, radius
HOUSE_EXTENT = (-.84, .84, -1.52, .66)  # House-local x/y extent: eaves, porch, path.
SOURCE_FOOTPRINT = .32                  # Source pedestal radius plus a margin.

def check_source_site(source):
    """Reject an independent source that would stand in a house, tree or shrub."""
    px, py = source["xy"]
    for village in VILLAGES:
        ox, oy = px-village["xy"][0], py-village["xy"][1]
        r = math.hypot(ox, oy)
        if r >= VILLAGE_RADIUS + SOURCE_FOOTPRINT:
            continue
        where = "Source " + source["name"] + " in village " + village["name"]
        if r > VILLAGE_RADIUS - SOURCE_FOOTPRINT:
            raise ValueError(where + " straddles the village boundary; move it inside or outside.")
        for hi, (hx, hy, rot, scale) in enumerate(HOUSE_LAYOUT):
            a = math.radians(rot)
            lx = ( math.cos(a)*(ox-hx) + math.sin(a)*(oy-hy)) / scale
            ly = (-math.sin(a)*(ox-hx) + math.cos(a)*(oy-hy)) / scale
            qx = max(HOUSE_EXTENT[0]-lx, 0, lx-HOUSE_EXTENT[1])
            qy = max(HOUSE_EXTENT[2]-ly, 0, ly-HOUSE_EXTENT[3])
            if math.hypot(qx, qy)*scale < SOURCE_FOOTPRINT:
                raise ValueError(where + " overlaps house " + str(hi+1) +
                                 " (its walls, porch or path); try another offset.")
        (tx, ty), tscale = VILLAGE_TREE
        if math.hypot(ox-tx, oy-ty) < .95*tscale + SOURCE_FOOTPRINT:
            raise ValueError(where + " is under the village tree; try another offset.")
        for (sx, sy), radius in VILLAGE_SHRUBS:
            if math.hypot(ox-sx, oy-sy) < radius + SOURCE_FOOTPRINT:
                raise ValueError(where + " overlaps a shrub; try another offset.")

def resolve_configuration(scenario):
    """Return explicit house assignments and additive source coordinates."""
    presets = {
        "untreated": ({}, []),
        "single_house": ({"B": [1]}, []),
        "two_houses": ({"B": [1, 3]}, []),
        "near_village": ({"B": "all"}, []),
        "far_village": ({"C": "all"}, []),
        "separate_source": ({}, [{"name": "S1", "village": SEPARATE_SOURCE_VILLAGE,
                                   "offset": SEPARATE_SOURCE_OFFSET,
                                   "strength": SEPARATE_SOURCE_STRENGTH}]),
        "manual": (MANUAL_TREATED_HOUSES, MANUAL_EXTERNAL_SOURCES),
    }
    if scenario not in presets:
        raise ValueError("Unknown SCENARIO. Choose one of: " + ", ".join(presets))
    if KERNEL not in ("exponential", "gaussian"):
        raise ValueError("KERNEL must be 'exponential' or 'gaussian'.")
    if SOURCE_UNIT not in ("houses", "village_centres"):
        raise ValueError("SOURCE_UNIT must be 'houses' or 'village_centres'.")
    unknown_styles = SPILLOVER_STYLES - {"shells", "terraces", "contours", "mist", "ground"}
    if unknown_styles or not SPILLOVER_STYLES:
        raise ValueError("SPILLOVER_STYLE: combine 'shells', 'terraces', 'contours', "
                         "'mist', 'ground' with '+'.")
    if (not all(math.isfinite(v) for v in (MIST_DENSITY,MIST_HEIGHT,MIST_GLOW,MIST_EDGE_FADE)) or
            MIST_DENSITY < 0 or MIST_HEIGHT <= .01 or MIST_GLOW < 0 or
            not 0 < MIST_EDGE_FADE < min(GROUND_SIZE)/2):
        raise ValueError("Mist density/glow must be nonnegative, height above .01, and edge fade within the board.")
    if (not FIELD_LEVELS or any(not v > 0 for v in FIELD_LEVELS) or
            list(FIELD_LEVELS) != sorted(FIELD_LEVELS, reverse=True) or
            len(set(FIELD_LEVELS)) != len(FIELD_LEVELS) or
            len(FIELD_LEVEL_COLOURS) < len(FIELD_LEVELS)):
        raise ValueError("FIELD_LEVELS must be distinct positive values, largest (innermost) "
                         "first, with one FIELD_LEVEL_COLOURS entry each.")
    if not (SHELL_HEIGHT > 0 and .01 <= SHELL_GRID <= .5 and 0 <= SHELL_OPACITY <= 1 and
            0 <= SHELL_RIM <= 1 and SHELL_GLOW >= 0 and CONTOUR_RADIUS > 0):
        raise ValueError("Check shell settings: positive height/radius, grid .01-.5, opacities in [0,1].")
    if BACKGROUND not in ("studio", "transparent"):
        raise ValueError("BACKGROUND must be 'studio' or 'transparent'.")
    if not (0 <= FRAME_MARGIN < .4 and SUN_STRENGTH >= 0 and SUN_ANGLE >= 0 and
            FILL_RATIO >= 0 and RIM_STRENGTH >= 0 and SKY_STRENGTH >= 0 and
            Vector(SUN_POSITION).z > 0):
        raise ValueError("Lights must be nonnegative, the sun above the horizon, margin in [0, .4).")
    if DECAY_LENGTH <= 0 or FIELD_REFERENCE <= 0 or not 0 <= MAX_FIELD_OPACITY <= 1:
        raise ValueError("Decay length/reference must be positive; opacity must be in [0,1].")
    positions = {v["name"]: v["xy"] for v in VILLAGES}
    if len(positions) != len(VILLAGES) or not HOUSE_LAYOUT:
        raise ValueError("Village names must be unique, and HOUSE_LAYOUT cannot be empty.")
    selected, external = presets[scenario]
    unknown = set(selected) - set(positions)
    if unknown:
        raise ValueError("Unknown village names: " + str(unknown))
    treated, sources = {}, []
    all_houses = set(range(1, len(HOUSE_LAYOUT)+1))
    for village, xy in positions.items():
        requested = selected.get(village, [])
        ids = all_houses.copy() if requested == "all" else set(requested)
        if not ids.issubset(all_houses) or any(type(i) is not int for i in ids):
            raise ValueError("House IDs must be integers from 1 to " + str(len(HOUSE_LAYOUT)))
        treated[village] = ids
        if SOURCE_UNIT == "houses":
            for house_id in sorted(ids):
                hx, hy, _, _ = HOUSE_LAYOUT[house_id-1]
                sources.append({"name": village + str(house_id),
                                "xy": (xy[0]+hx, xy[1]+hy),
                                "strength": HOUSE_SOURCE_STRENGTH,
                                "length": DECAY_LENGTH, "external": False})
        elif ids:
            if ids != all_houses:
                raise ValueError("Partial treatment needs SOURCE_UNIT = 'houses'.")
            sources.append({"name": village, "xy": xy,
                            "strength": VILLAGE_SOURCE_STRENGTH,
                            "length": DECAY_LENGTH, "external": False})
    for i, source in enumerate(external):
        if ("xy" in source) == ("village" in source):
            raise ValueError("An external source needs either xy or village + offset.")
        if "xy" in source:
            xy = tuple(source["xy"])
        else:
            if source["village"] not in positions:
                raise ValueError("Source village must be one of: " + ", ".join(positions))
            vx, vy = positions[source["village"]]
            dx, dy = source.get("offset", (0,0))
            xy = (vx+dx, vy+dy)
        sources.append({"name": source.get("name", "S"+str(i+1)), "xy": xy,
                        "strength": source.get("strength", 1.0),
                        "length": source.get("length", DECAY_LENGTH), "external": True})
    for source in sources:
        if (len(source["xy"]) != 2 or
                not all(math.isfinite(v) for v in (*source["xy"],source["strength"],source["length"])) or
                source["strength"] < 0 or source["length"] <= 0):
            raise ValueError("Sources need finite x,y, nonnegative strength and positive length.")
        if source["external"]:
            check_source_site(source)
    return treated, sources

def spillover_value(x, y, sources):
    """Same kernel as the ground shader, in ground-plane distance units."""
    total = 0.0
    for source in sources:
        d = math.hypot(x-source["xy"][0],y-source["xy"][1]) / source["length"]
        total += source["strength"] * math.exp(-d if KERNEL == "exponential" else -.5*d*d)
    return total


treated_houses, spillover_sources = resolve_configuration(SCENARIO)


# ----------------------- Scene and geometry tools ----------------------
scene = bpy.data.scenes.new("Village diorama | " + SCENARIO)
if bpy.context.window:
    bpy.context.window.scene = scene

def collection(name):
    result = bpy.data.collections.new(name)
    scene.collection.children.link(result)
    return result

land = collection("01 Landscape")
buildings = collection("02 Houses")
plants = collection("03 Trees and planting")
markers = collection("04 Village boundaries and markers")
studio = collection("05 Camera and lights")
field_collection = collection("06 Spillover field")

def colour(hex_colour):
    values = [int(hex_colour[i:i+2], 16) / 255 for i in (0, 2, 4)]
    return tuple(v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055)**2.4
                 for v in values) + (1.0,)

def material(name, light, dark=None, scale=5, roughness=.82,
             bump=.015, grain=None):
    """A two-scale procedural material: broad colour, fine surface texture."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    mat.diffuse_color = colour(light)
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    output = nodes.new("ShaderNodeOutputMaterial")
    shader = nodes.new("ShaderNodeBsdfPrincipled")
    shader.inputs["Base Color"].default_value = colour(light)
    shader.inputs["Roughness"].default_value = roughness
    links.new(shader.outputs["BSDF"], output.inputs["Surface"])
    if dark:
        coord = nodes.new("ShaderNodeTexCoord")
        vector = coord.outputs["Object"]
        if grain:
            stretch = nodes.new("ShaderNodeVectorMath")
            stretch.operation = "MULTIPLY"
            stretch.inputs[1].default_value = grain
            links.new(vector, stretch.inputs[0])
            vector = stretch.outputs["Vector"]
        noise = nodes.new("ShaderNodeTexNoise")
        noise.inputs["Scale"].default_value = scale
        noise.inputs["Detail"].default_value = 3.0
        links.new(vector, noise.inputs["Vector"])
        ramp = nodes.new("ShaderNodeValToRGB")
        ramp.color_ramp.elements[0].position = .18
        ramp.color_ramp.elements[0].color = colour(dark)
        ramp.color_ramp.elements[1].position = .82
        ramp.color_ramp.elements[1].color = colour(light)
        links.new(noise.outputs["Fac"], ramp.inputs["Fac"])
        links.new(ramp.outputs["Color"], shader.inputs["Base Color"])
        if bump:
            fine = nodes.new("ShaderNodeTexNoise")
            fine.inputs["Scale"].default_value = scale * 12
            fine.inputs["Detail"].default_value = 2
            links.new(vector, fine.inputs["Vector"])
            relief = nodes.new("ShaderNodeBump")
            relief.inputs["Strength"].default_value = .28
            relief.inputs["Distance"].default_value = bump
            links.new(fine.outputs["Fac"], relief.inputs["Height"])
            links.new(relief.outputs["Normal"], shader.inputs["Normal"])
    return mat

grass = material("Grass | mottled sage", "A4B77B", "758C53", 3, bump=.025)
soil = material("Soil | warm earth", "96704C", "634C39", 7, bump=.045)
topsoil = material("Topsoil | dark layer", "766043", "4E4230", 11, bump=.024)
yard = material("Village ground | sandy gravel", "CFC4A2", "A9A27A", 10, bump=.018)
stone = material("Stone | warm limestone", "BFB8A5", "959181", 12, bump=.025)
plasters = [material("Plaster " + str(i), hi, lo, 5, bump=.012)
            for i, (hi, lo) in enumerate([
                ("F3E7C9", "D4C7A8"), ("E9DEBF", "C9BC9A"),
                ("F3DBBE", "D6B79A")])]
wood = material("Timber | fine grain", "8A6846", "55402D", 4,
                bump=.009, grain=(7, 7, .5))
trim = material("Painted joinery | warm ivory", "F4E9CD", "DDD0AF", 8, bump=.006)
glass = material("Window panes | blue grey", "54767A", roughness=.24, bump=0)
shutter = material("Shutters | muted teal", "537D78", "365A57", 5, bump=.008)
roof_mats = [material("Clay tiles " + str(i), hi, lo, 8, bump=.012)
             for i, (hi, lo) in enumerate([
                 ("BF7454", "A96247"), ("CA8360", "B16D50"),
                 ("B96C4D", "9E5B43"), ("C57C59", "AE664A"),
                 ("D18A66", "B97556")])]
treated_roof_mats = [material("Treated roof tiles " + str(i), hi, lo, 8, bump=.012)
                    for i, (hi, lo) in enumerate([
                        ("398ECD", "2376AF"), ("459BD6", "3183BD"),
                        ("3186C4", "226FA7"), ("4096D0", "2A7AB4"),
                        ("50A1D9", "3788C0")])]
leaves = [material("Leaves " + str(i), hi, lo, 7, bump=.03)
          for i, (hi, lo) in enumerate([
              ("779951", "4E723C"), ("91A95F", "67884C"),
              ("A6B971", "7F9C58")])]
ring_mat = material("Village boundary | ivory", "F1E8CE", roughness=.9)
treated_mat = material("Treatment | blue", "2386C7", roughness=.65)
annotation_mat = material("Annotations | charcoal", "394A54", roughness=1)

def mesh_object(name, verts, faces, mat, coll, parent=None):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    coll.objects.link(obj)
    obj.parent = parent
    if mat:
        mesh.materials.append(mat)
    return obj

CUBE_VERTS = [(-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
              (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1)]
CUBE_FACES = [(0,3,2,1), (4,5,6,7), (0,1,5,4), (1,2,6,5),
              (2,3,7,6), (3,0,4,7)]

def box(name, loc, size, mat, coll, parent=None, bevel=.02):
    verts = [tuple(v[k] * size[k] / 2 for k in range(3)) for v in CUBE_VERTS]
    obj = mesh_object(name, verts, CUBE_FACES, mat, coll, parent)
    obj.location = loc
    if bevel:
        mod = obj.modifiers.new("Soft edges", "BEVEL")
        mod.width = bevel
        mod.segments = 2
    return obj

def cylinder(name, loc, radius, depth, mat, coll, parent=None, radius_top=None, n=32):
    rt = radius if radius_top is None else radius_top
    verts = [(r * math.cos(2*math.pi*i/n), r * math.sin(2*math.pi*i/n), z)
             for r, z in ((radius, -depth/2), (rt, depth/2)) for i in range(n)]
    faces = [tuple(reversed(range(n))), tuple(range(n, 2*n))]
    faces += [(i, (i+1)%n, (i+1)%n+n, i+n) for i in range(n)]
    obj = mesh_object(name, verts, faces, mat, coll, parent)
    obj.location = loc
    return obj

def stroke(name, coords, radius, mat, coll, parent=None, closed=False):
    data = bpy.data.curves.new(name, "CURVE")
    data.dimensions = "3D"
    data.bevel_depth = radius
    data.bevel_resolution = 2
    spline = data.splines.new("POLY")
    spline.points.add(len(coords)-1)
    for point, co in zip(spline.points, coords):
        point.co = (*co, 1)
    spline.use_cyclic_u = closed
    obj = bpy.data.objects.new(name, data)
    coll.objects.link(obj)
    obj.parent = parent
    data.materials.append(mat)
    return obj

def empty(name, loc, coll, parent=None):
    obj = bpy.data.objects.new(name, None)
    coll.objects.link(obj)
    obj.parent = parent
    obj.location = loc
    obj.empty_display_size = .15
    return obj

def pebble(name, loc, size, mat, coll, parent=None, smooth=True):
    mesh = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=2, radius=1)
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new(name, mesh)
    coll.objects.link(obj)
    obj.parent = parent
    obj.location, obj.scale = loc, size
    mesh.materials.append(mat)
    for p in mesh.polygons:
        p.use_smooth = smooth
    return obj


# ---------------------------- Landscape --------------------------------
gx, gy = GROUND_SIZE
box("Exposed soil block", (0,0,-.43), (gx,gy,.70), soil, land, bevel=.16)
box("Dark topsoil layer", (0,0,-.10), (gx-.025,gy-.025,.15), topsoil, land, bevel=.06)
box("Grass cap", (0,0,-.017), (gx-.05,gy-.05,.075), grass, land, bevel=.04)

rng = random.Random(SEED)
# Small embedded stones make the cut soil edges read as real material.
for edge in range(4):
    for i in range(18):
        z = rng.uniform(-.62,-.20)
        if edge < 2:
            loc = (rng.uniform(-gx/2+.3,gx/2-.3), (-1 if edge==0 else 1)*(gy/2-.015), z)
        else:
            loc = ((-1 if edge==2 else 1)*(gx/2-.015), rng.uniform(-gy/2+.3,gy/2-.3), z)
        r = rng.uniform(.035,.075)
        pebble("Stone in soil", loc, (r, r*.65, r*.48), stone, land)


# ------------------------------ Houses ---------------------------------
def window(parent, x, y=-.455, z=.60):
    box("Window recess", (x,y-.012,z), (.285,.025,.32), wood, buildings, parent, .008)
    box("Blue window panes", (x,y-.031,z), (.245,.016,.275), glass, buildings, parent, .004)
    for dx in (-.146,.146):
        box("Window frame", (x+dx,y-.047,z), (.031,.032,.35), trim, buildings, parent, .004)
    for dz in (-.16,.16):
        box("Window frame", (x,y-.047,z+dz), (.32,.032,.032), trim, buildings, parent, .004)
    box("Window mullion", (x,y-.054,z), (.020,.023,.29), trim, buildings, parent, .003)
    box("Window crossbar", (x,y-.054,z), (.26,.023,.020), trim, buildings, parent, .003)
    box("Projecting sill", (x,y-.085,z-.18), (.36,.15,.040), stone, buildings, parent, .008)
    for dx in (-.23,.23):
        box("Shutter", (x+dx,y-.025,z), (.12,.045,.32), shutter, buildings, parent, .008)
        for dz in (-.09,0,.09):
            box("Shutter slat", (x+dx,y-.052,z+dz), (.10,.018,.019), wood, buildings, parent, .002)

def tiled_roof(parent, local_rng, palette):
    # Two roof halves. Each contains individually modelled overlapping tiles,
    # batched into one mesh so the scene stays responsive.
    half_span, rise, ridge = .62, .48, 1.43
    angle = math.atan2(rise, half_span)
    slope_length = math.hypot(half_span,rise)
    for side in (-1,1):
        verts, faces, indices = [], [], []
        rotation = Matrix.Rotation(-side*angle, 3, "X")
        nx, nr = 11, 7
        for row in range(nr):
            t = (row+.5)/nr
            for col in range(nx):
                width = 1.61/nx
                centre = Vector((-0.805+(col+.5)*width,
                                 side*half_span*t, ridge-rise*t))
                size = (width-.007, slope_length/nr*1.09, .030)
                start = len(verts)
                verts.extend(tuple(centre+rotation@Vector(tuple(v[k]*size[k]/2 for k in range(3))))
                             for v in CUBE_VERTS)
                faces.extend(tuple(start+j for j in face) for face in CUBE_FACES)
                indices.extend([local_rng.randrange(len(palette))]*6)
        obj = mesh_object("Individual clay roof tiles", verts, faces, None, buildings, parent)
        for mat in palette:
            obj.data.materials.append(mat)
        for face, idx in zip(obj.data.polygons, indices):
            face.material_index = idx
        bevel = obj.modifiers.new("Tile edge rounding", "BEVEL")
        bevel.width, bevel.segments = .006, 2
        stroke("Eave fascia", [(-.83,side*.62,.93),(.83,side*.62,.93)], .035, wood, buildings, parent)
        for x in (-.805,.805):
            stroke("Gable fascia", [(x,0,1.43),(x,side*.62,.94)], .026, trim, buildings, parent)
    for i in range(10):
        tile = cylinder("Rounded ridge tile", (-.765+(i+.5)*.153,0,1.44),
                        .047,.161,palette[i%len(palette)],buildings,parent,n=12)
        tile.rotation_euler[1] = math.pi/2

def house(name, loc, rotation, scale, variant, parent, treated=False):
    h = empty(name, loc, buildings, parent)
    h.rotation_euler[2] = math.radians(rotation)
    h.scale = (scale,)*3
    h["treated"] = treated
    local_rng = random.Random(SEED+100+variant)
    palette = treated_roof_mats if treated else roof_mats
    wallmat = plasters[variant%len(plasters)]
    box("Stone foundation", (0,0,.12), (1.44,.98,.20), stone, buildings, h, .035)
    box("Textured plaster walls", (0,0,.61), (1.35,.90,.94), wallmat, buildings, h, .025)
    # Gable infill at both ends of the ridge.
    verts = [(x,y,z) for x in (-.675,.675) for y,z in ((-.45,1.08),(.45,1.08),(0,1.43))]
    mesh_object("Plaster gable ends",verts,[(0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)],wallmat,buildings,h)
    tiled_roof(h,local_rng,palette)
    # Entrance on the long front facade, with visible planking and a doorstep.
    box("Door surround", (0,-.473,.40), (.29,.08,.61), trim, buildings, h, .01)
    box("Timber door", (0,-.519,.39), (.23,.026,.55), wood, buildings, h, .006)
    for x in (-.07,0,.07):
        box("Door plank seam", (x,-.537,.39), (.007,.008,.51), soil, buildings, h, .001)
    pebble("Door handle",(.075,-.553,.40),(.018,.017,.018),stone,buildings,h)
    box("Upper doorstep",(0,-.58,.13),(.42,.30,.15),stone,buildings,h,.015)
    box("Lower doorstep",(0,-.75,.055),(.48,.20,.08),stone,buildings,h,.015)
    window(h,-.43)
    window(h,.43)
    # Matching window on the right end wall.
    end = empty("End-wall window",(.676,0,0),buildings,h)
    end.rotation_euler[2] = math.pi/2
    window(end,0,y=0,z=.60)
    # Small porch roof, with two slender wooden supports.
    awning = box("Porch canopy",(0,-.70,.89),(.53,.47,.055),palette[2],buildings,h,.012)
    awning.rotation_euler[0] = .12
    for x in (-.23,.23):
        box("Porch post",(x,-.88,.45),(.035,.035,.80),wood,buildings,h,.005)
    # Short chimney emerging through the roof.
    box("Plastered chimney",(-.42,.23,1.40),(.17,.18,.45),wallmat,buildings,h,.012)
    box("Chimney cap",(-.42,.23,1.64),(.22,.23,.045),stone,buildings,h,.008)
    cylinder("Chimney pot",(-.42,.23,1.70),.044,.10,roof_mats[0],buildings,h,n=12)
    for i in range(3):
        pebble("Path stepping stone",(0,-1.02-i*.19,.039),(.19,.115,.035),stone,buildings,h,smooth=False)
    if treated and SHOW_TREATMENT_FLAGS:
        stroke("House treatment flagpole",[(.43,0,1.43),(.43,0,1.97)],.016,wood,markers,h)
        mesh_object("House treatment pennant",[(.43,-.015,1.96),(.78,-.015,1.85),(.43,-.015,1.73)],
                    [(0,1,2)],treated_mat,markers,h)
    return h


# ------------------------- Trees and planting --------------------------
def tree(name, loc, size, seed, parent=None):
    tr = empty(name,loc,plants,parent)
    tr.scale = (size,)*3
    r = random.Random(seed)
    cylinder("Tapered trunk",(0,0,.63),.105,1.26,wood,plants,tr,radius_top=.054,n=10)
    for i in range(5):
        a = i*2*math.pi/5 + .3
        tip = (.37*math.cos(a),.37*math.sin(a),1.35+r.uniform(-.1,.12))
        stroke("Branch",[(0,0,.65),(tip[0]*.48,tip[1]*.48,1.02),tip],.038,wood,plants,tr)
        pebble("Leafy crown",(tip[0],tip[1],tip[2]+.24),
               (.49+r.random()*.10,.48+r.random()*.10,.49+r.random()*.12),
               leaves[i%3],plants,tr)
    pebble("Top crown",(.02,.015,1.90),(.48,.45,.46),leaves[2],plants,tr)
    for i in range(3):
        a = i*2.1
        stroke("Visible root",[(0,0,.13),(.20*math.cos(a),.20*math.sin(a),.026)],.033,wood,plants,tr)
    return tr

for vi, village in enumerate(VILLAGES):
    x,y = village["xy"]
    root = empty("Village " + village["name"],(x,y,0),markers)
    root["centre_x"],root["centre_y"] = x,y
    selected_houses = treated_houses[village["name"]]
    root["treated_house_ids"] = ",".join(str(i) for i in sorted(selected_houses))
    cylinder("Village " + village["name"] + " clearing",(0,0,.018),VILLAGE_RADIUS,.020,yard,land,root,n=128)
    if SHOW_VILLAGE_RINGS:
        points = [(VILLAGE_RADIUS*math.cos(t*math.pi/64),
                   VILLAGE_RADIUS*math.sin(t*math.pi/64),.050) for t in range(128)]
        focal = village.get("role") == "focal"
        stroke("Village " + village["name"] + " boundary",points,.030 if focal else .024,
               annotation_mat if focal else ring_mat,markers,root,True)
    for hi,(hx,hy,rot,scale) in enumerate(HOUSE_LAYOUT):
        h = house("House " + village["name"] + str(hi+1),(hx,hy,.035),rot,scale,
                  vi*len(HOUSE_LAYOUT)+hi,root,treated=(hi+1 in selected_houses))
        h["spillover_S"] = spillover_value(x+hx,y+hy,spillover_sources)
    (tx,ty),tscale = VILLAGE_TREE
    tree("Village tree " + village["name"],(tx,ty,.03),tscale,SEED+vi*11,root)
    # Low shrubs along the edge, leaving the foreground entrance clear.
    for j,((sx,sy),radius) in enumerate(VILLAGE_SHRUBS):
        pebble("Village shrub",(sx,sy,.15),(radius,radius*.83,radius*.83),leaves[j%3],plants,root)

# Sparse scenery around the perimeter preserves space for later distance lines.
tree("Landscape tree west",(-5.7,-2.7,.025),1.03,SEED+200)
tree("Landscape tree east",(5.7,4.0,.025),1.08,SEED+201)
for i in range(70):
    x,y = rng.uniform(-gx/2+.3,gx/2-.3),rng.uniform(-gy/2+.3,gy/2-.3)
    if any(math.hypot(x-v["xy"][0],y-v["xy"][1]) < VILLAGE_RADIUS+.25 for v in VILLAGES):
        continue
    if i%4 == 0:
        pebble("Landscape stone",(x,y,.05),(.13,.09,.07),stone,land,smooth=False)
    else:
        for j in range(3):
            angle = rng.random()*2*math.pi
            height = rng.uniform(.07,.15)
            verts = [(x-.015,y,.025),(x+.015,y,.025),
                     (x+.06*math.cos(angle),y+.06*math.sin(angle),height)]
            mesh_object("Grass blade",verts,[(0,1,2)],leaves[2],plants)


# ---------------------- Sources and spillover field --------------------
for source in spillover_sources:
    if not source["external"]:
        continue
    x,y = source["xy"]
    root = empty("Independent source " + source["name"],(x,y,.035),markers)
    root["strength"] = source["strength"]
    root["decay_length"] = source["length"]
    cylinder("Source pedestal",(0,0,.065),.23,.13,stone,markers,root,n=32)
    cylinder("Source post",(0,0,.53),.065,.86,wood,markers,root,n=16)
    pebble("Blue source marker",(0,0,1.02),(.17,.17,.17),treated_mat,markers,root)
    mesh_object("Source pennant",[(.02,0,.88),(.48,0,.76),(.02,0,.62)],
                [(0,1,2)],treated_mat,markers,root)
    points = [(.32*math.cos(t*math.pi/24),.32*math.sin(t*math.pi/24),.025) for t in range(48)]
    stroke("Source ground marker",points,.025,treated_mat,markers,root,True)

def kernel_shader_nodes(nodes, links, sources):
    """Shared world-coordinate kernel calculation for ground and mist shaders."""
    def calculate(operation, a, b=None):
        node = nodes.new("ShaderNodeMath")
        node.operation = operation
        for i,value in enumerate((a,) if b is None else (a,b)):
            if isinstance(value,(int,float)):
                node.inputs[i].default_value = value
            else:
                links.new(value,node.inputs[i])
        return node.outputs[0]

    geometry = nodes.new("ShaderNodeNewGeometry")
    split = nodes.new("ShaderNodeSeparateXYZ")
    planar = nodes.new("ShaderNodeCombineXYZ")
    links.new(geometry.outputs["Position"],split.inputs["Vector"])
    links.new(split.outputs["X"],planar.inputs["X"])
    links.new(split.outputs["Y"],planar.inputs["Y"])
    planar.inputs["Z"].default_value = 0
    total = 0.0
    for source in sources:
        distance = nodes.new("ShaderNodeVectorMath")
        distance.operation = "DISTANCE"
        links.new(planar.outputs["Vector"],distance.inputs[0])
        distance.inputs[1].default_value = (*source["xy"],0)
        scaled = calculate("DIVIDE",distance.outputs["Value"],source["length"])
        exponent = (calculate("MULTIPLY",scaled,-1.0) if KERNEL == "exponential" else
                    calculate("MULTIPLY",calculate("MULTIPLY",scaled,scaled),-.5))
        term = calculate("MULTIPLY",calculate("POWER",math.e,exponent),source["strength"])
        total = calculate("ADD",total,term)
    return total, split, calculate

def spillover_material(sources):
    """Exact planar kernels, summed BEFORE applying the fixed colour mapping.

    Camera-only colour overlay: no inverse-square lamp falloff, cast shadows,
    coloured bounce light or invented barriers to spatial spillover.
    """
    mat = bpy.data.materials.new("Spillover | additive " + KERNEL)
    mat.use_nodes = True
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    total, position, calculate = kernel_shader_nodes(nodes,links,sources)
    exponent = calculate("MULTIPLY",total,-1.0/FIELD_REFERENCE)
    opacity = calculate("MULTIPLY",calculate("SUBTRACT",1.0,calculate("POWER",math.e,exponent)),
                        MAX_FIELD_OPACITY)
    rays = nodes.new("ShaderNodeLightPath")
    opacity = calculate("MULTIPLY",opacity,rays.outputs["Is Camera Ray"])
    transparent = nodes.new("ShaderNodeBsdfTransparent")
    emission = nodes.new("ShaderNodeEmission")
    emission.inputs["Color"].default_value = colour(FIELD_COLOUR)
    emission.inputs["Strength"].default_value = 1.0
    mix = nodes.new("ShaderNodeMixShader")
    links.new(opacity,mix.inputs[0])
    links.new(transparent.outputs[0],mix.inputs[1])
    links.new(emission.outputs[0],mix.inputs[2])
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(mix.outputs[0],output.inputs["Surface"])
    return mat

def mist_material(sources):
    """A translucent blue volume driven by the SAME additive horizontal field.

    There is one containing mesh and one density sum, not separate overlapping
    gas objects. Density and glow are zero for non-camera rays so the decorative
    mist does not cast shadows or add blue bounce light to untreated roofs.
    Appearance is still view-dependent because volume opacity integrates along
    the ray. This is an illustration of influence, not a quantitative legend.
    """
    mat = bpy.data.materials.new("Spillover mist | additive " + KERNEL)
    mat.use_nodes = True
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    total, position, calculate = kernel_shader_nodes(nodes,links,sources)

    # Flat at ground level; smooth fade around and above roof height.
    height_ratio = calculate("DIVIDE",position.outputs["Z"],MIST_HEIGHT)
    height_power = calculate("POWER",height_ratio,4.0)
    vertical = calculate("POWER",math.e,calculate("MULTIPLY",height_power,-1.0))

    # The mathematical field has an infinite tail. A smooth display crop keeps
    # that tail from ending at a visible rectangular wall at the mesh boundary.
    # This crop is independent of source positions and is fixed across panels.
    edge_masks = []
    for axis,extent in (("X",GROUND_SIZE[0]/2),("Y",GROUND_SIZE[1]/2)):
        clearance = calculate("SUBTRACT",extent,calculate("ABSOLUTE",position.outputs[axis]))
        t = calculate("MINIMUM",1.0,calculate("MAXIMUM",0.0,
                       calculate("DIVIDE",clearance,MIST_EDGE_FADE)))
        smooth = calculate("MULTIPLY",calculate("MULTIPLY",t,t),
                           calculate("SUBTRACT",3.0,calculate("MULTIPLY",2.0,t)))
        edge_masks.append(smooth)
    crop = calculate("MULTIPLY",edge_masks[0],edge_masks[1])
    density = calculate("MULTIPLY",calculate("MULTIPLY",total,vertical),MIST_DENSITY)
    density = calculate("MULTIPLY",density,crop)
    rays = nodes.new("ShaderNodeLightPath")
    density = calculate("MULTIPLY",density,rays.outputs["Is Camera Ray"])

    volume = nodes.new("ShaderNodeVolumePrincipled")
    volume.inputs["Color"].default_value = colour(MIST_COLOUR)
    volume.inputs["Anisotropy"].default_value = 0.0
    volume.inputs["Emission Color"].default_value = colour(MIST_COLOUR)
    volume.inputs["Blackbody Intensity"].default_value = 0.0
    # No simulation grids are used; all density comes from the procedural nodes.
    density_attribute = volume.inputs.get("Density Attribute")
    if density_attribute is not None:
        density_attribute.default_value = ""
    colour_attribute = volume.inputs.get("Color Attribute")
    if colour_attribute is not None:
        colour_attribute.default_value = ""
    links.new(density,volume.inputs["Density"])
    links.new(calculate("MULTIPLY",density,MIST_GLOW),volume.inputs["Emission Strength"])
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(volume.outputs["Volume"],output.inputs["Volume"])
    return mat

FIELD_INSET = .08    # Keep shells and contours on the flat top of the grass cap.

def field_grid():
    """Spillover S(x,y) on a regular grid covering the top of the board.

    No edge crop: levels that reach the board edge are simply cut there, as if
    the landscape continued beyond the diorama, so contour positions are exact.
    """
    wx, wy = GROUND_SIZE[0]/2 - FIELD_INSET, GROUND_SIZE[1]/2 - FIELD_INSET
    nx = int(math.ceil(2*wx/SHELL_GRID)) + 1
    ny = int(math.ceil(2*wy/SHELL_GRID)) + 1
    xs = [-wx + 2*wx*i/(nx-1) for i in range(nx)]
    ys = [-wy + 2*wy*j/(ny-1) for j in range(ny)]
    values = [[spillover_value(x,y,spillover_sources) for x in xs] for y in ys]
    return xs, ys, values

def level_geometry(level, xs, ys, values, base):
    """Marching squares: the region S > level as a dome mesh, plus its outline.

    Inside grid points are lifted to base + SHELL_HEIGHT*sqrt(ln(S/level)), the
    iso-surface of S(x,y)*exp(-(z/SHELL_HEIGHT)^2) = level. Contour crossings are
    interpolated along grid edges and lie at z = base, so each dome meets the
    ground exactly along the planar contour S(x,y) = level.
    Returns verts, faces and contour chains as (points, closed) pairs; chains
    that leave the board are open and end at the board edge.
    """
    nx, ny = len(xs), len(ys)
    def border(i, j):
        return {tag for tag,hit in (("x0",i==0),("x1",i==nx-1),("y0",j==0),("y1",j==ny-1)) if hit}
    index, verts, faces, on_border = {}, [], [], []
    def grid_point(i, j):
        if (i,j) not in index:
            index[(i,j)] = len(verts)
            verts.append((xs[i], ys[j],
                          base + SHELL_HEIGHT*math.sqrt(math.log(values[j][i]/level))))
            on_border.append(border(i, j))
        return index[(i,j)]
    def crossing(a, b):
        key = (min(a,b), max(a,b))
        if key not in index:
            fa, fb = values[a[1]][a[0]]-level, values[b[1]][b[0]]-level
            t = fa/(fa-fb)
            index[key] = len(verts)
            verts.append((xs[a[0]] + t*(xs[b[0]]-xs[a[0]]),
                          ys[a[1]] + t*(ys[b[1]]-ys[a[1]]), base))
            on_border.append(border(*a) & border(*b))
        return index[key]
    for j in range(ny-1):
        for i in range(nx-1):
            corners = ((i,j), (i+1,j), (i+1,j+1), (i,j+1))
            inside = [values[c[1]][c[0]] > level for c in corners]
            if not any(inside):
                continue
            polygon = []
            for k in range(4):
                if inside[k]:
                    polygon.append(grid_point(*corners[k]))
                if inside[k] != inside[(k+1)%4]:
                    polygon.append(crossing(corners[k], corners[(k+1)%4]))
            faces.append(tuple(polygon))
    # Mesh edges used by exactly one face are the contour, except those lying
    # along the board edge, which are where the dome is cut off.
    uses = {}
    for face in faces:
        for k in range(len(face)):
            edge = (min(face[k-1],face[k]), max(face[k-1],face[k]))
            uses[edge] = uses.get(edge, 0) + 1
    neighbours = {}
    for (a,b),count in uses.items():
        if count == 1 and not (on_border[a] & on_border[b]):
            neighbours.setdefault(a, []).append(b)
            neighbours.setdefault(b, []).append(a)
    chains, seen = [], set()
    ends = [v for v,n in neighbours.items() if len(n) == 1]
    for start in ends + list(neighbours):     # Open chains first, then loops.
        if start in seen:
            continue
        chain, current = [start], start
        seen.add(start)
        while True:
            options = [v for v in neighbours[current] if v not in seen]
            if not options:
                break
            current = options[0]
            seen.add(current)
            chain.append(current)
        if len(chain) > 2:
            closed = start not in ends and chain[0] in neighbours[chain[-1]]
            chains.append(([verts[v] for v in chain], closed))
    return verts, faces, chains

def translucent_material(name, hex_colour, kind, top=1.0):
    """Camera-only translucent blue for shells ("shell"), terrace walls ("wall")
    and terrace floors ("floor").

    As with the mist, non-camera rays see pure transparency, so these surfaces
    cast no shadows and add no blue bounce light to roofs.
    """
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    mat.diffuse_color = (*colour(hex_colour)[:3], .20)   # Translucent in viewport.
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    def calculate(operation, a, b):
        node = nodes.new("ShaderNodeMath")
        node.operation = operation
        node.use_clamp = True
        for i,value in enumerate((a,b)):
            if isinstance(value,(int,float)):
                node.inputs[i].default_value = value
            else:
                links.new(value,node.inputs[i])
        return node.outputs[0]
    rays = nodes.new("ShaderNodeLightPath")
    camera_only = rays.outputs["Is Camera Ray"]
    if kind == "shell":
        # Faint face-on, stronger at grazing angles: the silhouette of each
        # dome reads as a surface rather than as a flat tint.
        weight = nodes.new("ShaderNodeLayerWeight")
        weight.inputs["Blend"].default_value = .45
        opacity = calculate("ADD", SHELL_OPACITY,
                            calculate("MULTIPLY", calculate("POWER", weight.outputs["Facing"], 2.0),
                                      SHELL_RIM))
    elif kind == "wall":
        # Densest at the foot, thinning upward but still visible at the top,
        # so the step height of each terrace can be read.
        geometry = nodes.new("ShaderNodeNewGeometry")
        height = nodes.new("ShaderNodeSeparateXYZ")
        links.new(geometry.outputs["Position"], height.inputs["Vector"])
        fade = calculate("SUBTRACT", 1.0, calculate("DIVIDE", height.outputs["Z"], top))
        opacity = calculate("MULTIPLY", calculate("ADD", .35,
                            calculate("MULTIPLY", calculate("POWER", fade, 1.3), .65)),
                            TERRACE_OPACITY)
    else:
        opacity = TERRACE_FLOOR
    opacity = calculate("MULTIPLY", opacity, camera_only)
    surface = nodes.new("ShaderNodeBsdfPrincipled")
    # Mostly self-lit: a darker diffuse base means the blue comes from the fixed
    # emission, so the field colours barely change with the studio lighting.
    surface.inputs["Base Color"].default_value = tuple(c*.45 for c in colour(hex_colour)[:3]) + (1,)
    surface.inputs["Roughness"].default_value = .28
    # Weak reflections: stronger lighting would otherwise put a white sheen
    # over the villages inside the domes.
    surface.inputs["Specular IOR Level"].default_value = .12
    surface.inputs["Emission Color"].default_value = colour(hex_colour)
    surface.inputs["Emission Strength"].default_value = SHELL_GLOW
    transparent = nodes.new("ShaderNodeBsdfTransparent")
    mix = nodes.new("ShaderNodeMixShader")
    links.new(opacity, mix.inputs[0])
    links.new(transparent.outputs[0], mix.inputs[1])
    links.new(surface.outputs[0], mix.inputs[2])
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(mix.outputs[0], output.inputs["Surface"])
    return mat

if SHOW_SPILLOVER and spillover_sources and SPILLOVER_STYLES & {"shells","terraces","contours"}:
    xs, ys, grid_values = field_grid()
    for k,level in enumerate(FIELD_LEVELS):
        # Tiny lift per level avoids coincident faces where contours nearly touch.
        base = .040 + .002*(len(FIELD_LEVELS)-k)
        verts, faces, chains = level_geometry(level, xs, ys, grid_values, base)
        if not faces:
            continue          # This level is not reached anywhere in this panel.
        tag = "S >= " + format(level, "g")
        hue = FIELD_LEVEL_COLOURS[k]
        if "shells" in SPILLOVER_STYLES:
            dome = mesh_object("Spillover shell | " + tag, verts, faces,
                               translucent_material("Spillover shell | " + tag, hue, "shell"),
                               field_collection)
            dome.data.polygons.foreach_set("use_smooth", [True]*len(dome.data.polygons))
            dome["field_level"] = level
            dome["description"] = ("Iso-surface of S(x,y)*exp(-(z/H)^2); footprint is the "
                                   "ground contour S(x,y) = level.")
        if "terraces" in SPILLOVER_STYLES:
            # Stepped glass: a wall on each contour, taller for inner levels,
            # over a floor tint that accumulates towards the sources.
            top = base + TERRACE_STEP*(len(FIELD_LEVELS)-k)
            floor = mesh_object("Spillover terrace floor | " + tag,
                                [(x,y,base) for x,y,_ in verts], faces,
                                translucent_material("Spillover terrace floor | " + tag,
                                                     hue, "floor"), field_collection)
            floor["field_level"] = level
            wall_mat = translucent_material("Spillover terrace wall | " + tag, hue, "wall", top)
            for n,(chain,closed) in enumerate(chains):
                m = len(chain)
                wall = mesh_object("Spillover terrace wall | " + tag + " | " + str(n+1),
                                   [(x,y,base) for x,y,_ in chain] + [(x,y,top) for x,y,_ in chain],
                                   [(i,(i+1)%m,(i+1)%m+m,i+m) for i in range(m if closed else m-1)],
                                   wall_mat, field_collection)
                wall["field_level"] = level
        if "contours" in SPILLOVER_STYLES:
            line = material("Spillover contour | " + tag, hue, roughness=.6)
            for n,(chain,closed) in enumerate(chains):
                ring = stroke("Spillover contour | " + tag + " | " + str(n+1),
                              [(x,y,.046) for x,y,_ in chain], CONTOUR_RADIUS, line,
                              field_collection, closed=closed)
                ring["field_level"] = level

if SHOW_SPILLOVER and spillover_sources and "ground" in SPILLOVER_STYLES:
    # One overlay avoids alpha-stacking separately normalised source disks.
    # It lies above grass and village clearings, below the houses and markers.
    half_x,half_y = gx/2-.10,gy/2-.10
    overlay = mesh_object("Spillover gradient | view with F12",
                          [(-half_x,-half_y,.034),(half_x,-half_y,.034),
                           (half_x,half_y,.034),(-half_x,half_y,.034)],
                          [(0,1,2,3)],spillover_material(spillover_sources),field_collection)
    overlay.display_type = "WIRE"  # Avoid an opaque rectangle in solid viewport mode.

if SHOW_SPILLOVER and spillover_sources and "mist" in SPILLOVER_STYLES:
    # At 2.5 heights the vertical envelope is exp(-2.5^4), effectively zero.
    bottom,top = .025,2.5*MIST_HEIGHT
    cloud = box("Spillover mist | render with F12",(0,0,(bottom+top)/2),
                (gx,gy,top-bottom),mist_material(spillover_sources),field_collection,bevel=0)
    cloud.display_type = "WIRE"
    cloud["density_scale"] = MIST_DENSITY
    cloud["vertical_height"] = MIST_HEIGHT
    cloud["description"] = "Additive planar field with a decorative vertical fade and fixed edge crop."


# --------------------------- Camera and light --------------------------
def backdrop_material():
    """Cool grey studio sweep that darkens towards the frame corners."""
    mat = bpy.data.materials.new("Backdrop | cool grey vignette")
    mat.use_nodes = True
    mat.diffuse_color = colour(BACKDROP_COLOURS[0])
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    geometry = nodes.new("ShaderNodeNewGeometry")
    planar = nodes.new("ShaderNodeVectorMath")
    planar.operation = "MULTIPLY"
    planar.inputs[1].default_value = (1, 1, 0)
    links.new(geometry.outputs["Position"], planar.inputs[0])
    radius = nodes.new("ShaderNodeVectorMath")
    radius.operation = "LENGTH"
    links.new(planar.outputs["Vector"], radius.inputs[0])
    ramp = nodes.new("ShaderNodeValToRGB")
    reach = math.hypot(*GROUND_SIZE) / 2          # Board corner distance.
    ramp.color_ramp.interpolation = "EASE"
    ramp.color_ramp.elements[0].position = 0.0
    ramp.color_ramp.elements[0].color = colour(BACKDROP_COLOURS[0])
    ramp.color_ramp.elements[1].position = 1.0
    ramp.color_ramp.elements[1].color = colour(BACKDROP_COLOURS[1])
    # Map radius so the fade starts just outside the board and ends off-frame.
    scale = nodes.new("ShaderNodeMapRange")
    scale.inputs["From Min"].default_value = reach * .75
    scale.inputs["From Max"].default_value = reach * 2.1
    links.new(radius.outputs["Value"], scale.inputs["Value"])
    links.new(scale.outputs["Result"], ramp.inputs["Fac"])
    shader = nodes.new("ShaderNodeBsdfPrincipled")
    shader.inputs["Roughness"].default_value = 1.0
    links.new(ramp.outputs["Color"], shader.inputs["Base Color"])
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(shader.outputs["BSDF"], output.inputs["Surface"])
    return mat

floor = box("Studio floor",(0,0,-.83),(200,200,.10),backdrop_material(),studio,bevel=0)
floor.is_shadow_catcher = BACKGROUND == "transparent"

sun_vector = Vector(SUN_POSITION).normalized()

# Skylight: Blender's built-in Sky Texture (no image files), with its sun disc
# switched off because the Sun lamp below provides direct sunlight. The sky's
# sun is placed in the same direction so the soft sky gradient agrees with it.
world = bpy.data.worlds.new("Diorama daylight | sky texture")
world.use_nodes = True
sky = world.node_tree.nodes.new("ShaderNodeTexSky")
available = [item.identifier for item in sky.bl_rna.properties["sky_type"].enum_items]
for sky_type in ("MULTIPLE_SCATTERING", "NISHITA", "SINGLE_SCATTERING", "HOSEK_WILKIE"):
    if sky_type in available:
        sky.sky_type = sky_type
        break
if hasattr(sky, "sun_disc"):
    sky.sun_disc = False
if hasattr(sky, "sun_elevation"):
    sky.sun_elevation = math.asin(sun_vector.z)
    sky.sun_rotation = math.atan2(-sun_vector.x, sun_vector.y) % (2*math.pi)
background = world.node_tree.nodes["Background"]
background.inputs["Strength"].default_value = SKY_STRENGTH
world.node_tree.links.new(sky.outputs["Color"], background.inputs["Color"])
scene.world = world

def aim(obj, target=Vector((0,0,0))):
    obj.rotation_euler = (target-obj.location).to_track_quat('-Z','Y').to_euler()

def area_light(name, loc, energy, size, tint, shadows=True):
    data = bpy.data.lights.new(name,"AREA")
    data.energy, data.shape, data.size = energy,"DISK",size
    data.color = colour(tint)[:3]
    data.use_shadow = shadows
    obj = bpy.data.objects.new(name,data)
    studio.objects.link(obj)
    obj.location = loc
    aim(obj)
    return obj

sun_data = bpy.data.lights.new("Warm key sun","SUN")
sun_data.energy = SUN_STRENGTH
sun_data.angle = math.radians(SUN_ANGLE)
sun_data.color = colour("FFE8CC")[:3]
sun = bpy.data.objects.new("Warm key sun",sun_data)
studio.objects.link(sun)
sun.location = sun_vector * 15
aim(sun)
# Fill from the camera's right; large and soft so it lifts shadows without
# adding a second set of shadow edges. Energy is scaled from the key.
area_light("Cool fill",(7,1,9),SUN_STRENGTH*FILL_RATIO*260,9,"DCE8FF")
# Rim from behind the board, opposite the camera, low enough to graze roofs.
area_light("Cool rim",(-10,11,7),RIM_STRENGTH*650,5,"D2E2FF")

camera_data = bpy.data.cameras.new("Isometric camera")
camera_data.type = "ORTHO"
camera_data.lens = 50
camera_data.clip_end = 500
camera = bpy.data.objects.new("Isometric camera",camera_data)
studio.objects.link(camera)
target = Vector((0,0,.15))
camera.location = target + Vector((18,-18,18))
aim(camera, target)
scene.camera = camera

scene.render.resolution_x,scene.render.resolution_y = RESOLUTION
scene.render.resolution_percentage = 100

def frame_board():
    """Zoom and centre the orthographic camera on the board.

    Uses only the board, houses and planting, never flags, sources or the
    spillover display, so every scenario gets the identical frame.
    """
    scene.view_layers[0].update()                 # Resolve parented positions.
    rotation = camera.rotation_euler.to_matrix()
    right, up = rotation @ Vector((1,0,0)), rotation @ Vector((0,1,0))
    us, vs = [], []
    for coll in (land, buildings, plants):
        for obj in coll.all_objects:
            if obj.type not in ("MESH","CURVE"):
                continue
            for corner in obj.bound_box:
                point = obj.matrix_world @ Vector(corner)
                us.append(point.dot(right))
                vs.append(point.dot(up))
    width, height = max(us)-min(us), max(vs)-min(vs)
    x, y = RESOLUTION
    fit = max(width, height*x/y) if x >= y else max(width*y/x, height)
    camera_data.ortho_scale = fit / (1 - 2*FRAME_MARGIN)
    centre_u, centre_v = (max(us)+min(us))/2, (max(vs)+min(vs))/2
    camera.location += right*(centre_u - camera.location.dot(right))
    camera.location += up*(centre_v - camera.location.dot(up))

frame_board()

scene.render.engine = "CYCLES"
scene.cycles.device = "CPU"          # Works without configuring a GPU.
scene.cycles.samples = SAMPLES
scene.cycles.use_denoising = True
scene.cycles.seed = SEED
# A camera ray can cross the front and back of every nested shell, plus the
# ground overlay, before reaching a roof; the default limit would turn it black.
scene.cycles.transparent_max_bounces = max(scene.cycles.transparent_max_bounces,
                                           4*len(FIELD_LEVELS) + 8)
if SHOW_SPILLOVER and "mist" in SPILLOVER_STYLES:
    # Avoid relying on GPU settings; the previous CPU rendering workflow remains.
    # A smaller volume step resolves the smooth density around individual houses.
    if hasattr(scene.cycles,"volume_step_rate"):
        scene.cycles.volume_step_rate = .5
scene.render.image_settings.file_format = "PNG"
scene.render.image_settings.color_mode = "RGBA"
scene.render.film_transparent = BACKGROUND == "transparent"
scene.view_settings.view_transform = "AgX"
try:
    scene.view_settings.look = COLOUR_LOOK
except TypeError:
    print("Colour look not available here; using AgX default:", COLOUR_LOOK)
scene.view_settings.exposure = EXPOSURE

# Camera-facing labels. These optional annotations stay fixed between treatment
# configurations. Distances are computed from the same coordinates
# that place the villages; they are not read off the projected image.
def label(body, loc, size=.30):
    data = bpy.data.curves.new(body,"FONT")
    data.body = body
    data.size = size
    data.align_x = "CENTER"
    data.align_y = "CENTER"
    data.space_character = 1.05
    obj = bpy.data.objects.new(body,data)
    markers.objects.link(obj)
    obj.location = loc
    obj.rotation_euler = camera.rotation_euler.copy()
    data.materials.append(annotation_mat)
    return obj

if SHOW_LABELS:
    # A source standing inside a village is named in that village's label, which
    # sits in clear ground in front of the circle; any label near the marker
    # itself would overlap the houses around it.
    hosted = {v["name"]: [] for v in VILLAGES}
    for source in spillover_sources:
        if not source["external"]:
            continue
        x,y = source["xy"]
        host = [v["name"] for v in VILLAGES
                if math.hypot(x-v["xy"][0],y-v["xy"][1]) < VILLAGE_RADIUS]
        if host:
            hosted[host[0]].append(source["name"])
        else:
            label(source["name"] + " | source",(x,y,1.45),.25)
    for village in VILLAGES:
        x,y = village["xy"]
        text = village["name"] + " | " + village.get("role","")
        if hosted[village["name"]]:
            text += " + source " + ", ".join(hosted[village["name"]])
        # Offset toward the camera on the ground, outside the village circle.
        label(text,(x+1.80,y-1.80,.16),.32)

if SHOW_DISTANCE_GUIDES:
    by_name = {v["name"]: v for v in VILLAGES}
    if all(name in by_name for name in ("A","B","C")):
        a = Vector((*by_name["A"]["xy"],.069))
        for other,offset in (("B",1.30),("C",.50)):
            b = Vector((*by_name[other]["xy"],.069))
            delta = b-a
            distance = delta.length
            if distance == 0:
                continue
            direction = delta/distance
            # Start and end at village centres; buildings naturally occlude
            # hidden portions. The caption explicitly states centre distance.
            cursor = 0.0
            while cursor < distance:
                p = a+direction*cursor
                q = a+direction*min(cursor+.13,distance)
                stroke("Centre distance A-"+other,[tuple(p),tuple(q)],.012,annotation_mat,markers)
                cursor += .25
            midpoint = (a+b)/2
            normal = Vector((-direction.y,direction.x,0))
            text_position = midpoint+normal*offset
            text_position.z = .20
            if SHOW_LABELS:
                label("A-"+other+": "+format(distance,".2f"),tuple(text_position),.27)

# Pixel anchors for adding text in a document: x right, y down from the top-left
# corner of the rendered PNG, at RESOLUTION. Identical across scenarios except
# for sources, because the frame never depends on treatment.
def pixel(point):
    ndc = world_to_camera_view(scene, camera, Vector(point))
    return [round(ndc.x*RESOLUTION[0], 1), round((1-ndc.y)*RESOLUTION[1], 1)]

scene.view_layers[0].update()
anchors = {}
for village in VILLAGES:
    x,y = village["xy"]
    anchors["village " + village["name"] + " centre"] = pixel((x,y,.05))
    anchors["village " + village["name"] + " label spot"] = pixel((x+1.80,y-1.80,.16))
    for hi,(hx,hy,_,scale) in enumerate(HOUSE_LAYOUT):
        anchors["house " + village["name"] + str(hi+1) + " roof"] = pixel((x+hx,y+hy,.035+1.45*scale))
for source in spillover_sources:
    if source["external"]:
        anchors["source " + source["name"] + " marker"] = pixel((*source["xy"],1.2))
by_name = {v["name"]: v["xy"] for v in VILLAGES}
for other in ("B","C"):
    if "A" in by_name and other in by_name:
        mid = [(by_name["A"][k]+by_name[other][k])/2 for k in (0,1)]
        anchors["distance A-" + other + " midpoint"] = pixel((*mid,.07))
anchor_record = {"image_size": list(RESOLUTION), "origin": "top-left, y down",
                 "scenario": SCENARIO, "points": anchors}
scene["label_anchors"] = json.dumps(anchor_record)
if LABEL_ANCHORS_FILE:
    if LABEL_ANCHORS_FILE.startswith("//") and not bpy.data.filepath:
        print("Save the .blend first to write", LABEL_ANCHORS_FILE, "- anchors printed below instead.")
    else:
        with open(bpy.path.abspath(LABEL_ANCHORS_FILE), "w") as handle:
            json.dump(anchor_record, handle, indent=1)
        print("Label anchors written to", bpy.path.abspath(LABEL_ANCHORS_FILE))
print("Label anchors (pixels at", RESOLUTION, "from top-left):")
for key, value in anchors.items():
    print("  ", key.ljust(26), value)

# Make the initial viewport useful without waiting for a ray-traced preview.
# The full procedural textures and lighting appear in the F12 render.
for screen in bpy.data.screens:
    for area in screen.areas:
        if area.type == "VIEW_3D":
            space = area.spaces.active
            space.shading.type = "SOLID"
            space.shading.color_type = "MATERIAL"
            space.overlay.show_overlays = False
            space.region_3d.view_perspective = "CAMERA"
            space.region_3d.view_camera_zoom = 0

scene["description"] = "Blue roofs/flags = assigned treatment; blue shells/contours/mist = illustrative spillover."
scene["field_levels"] = list(FIELD_LEVELS)
scene["shell_mapping"] = "shell k = iso-surface S(x,y)*exp(-(z/height)^2) = level_k; footprint S(x,y) = level_k"
scene["shell_height"] = SHELL_HEIGHT
scene["seed"] = SEED
scene["scenario"] = SCENARIO
scene["source_unit"] = SOURCE_UNIT
scene["kernel"] = KERNEL
scene["decay_length"] = DECAY_LENGTH
scene["total_source_strength"] = sum(s["strength"] for s in spillover_sources)
scene["colour_mapping"] = "opacity = max_opacity * (1 - exp(-sum_of_kernels/reference))"
scene["spillover_style"] = SPILLOVER_STYLE
scene["mist_mapping"] = "density = density_scale * sum_of_kernels * exp(-(z/height)^4) * fixed_edge_crop"
scene["mist_density"] = MIST_DENSITY
scene["mist_height"] = MIST_HEIGHT
scene["distance_units"] = "Arbitrary scene units; guides connect village centres."
print("Village diorama ready. Use Render > Render Image (F12), then File > Save As.")
print("Scene:",scene.name,"| Objects:",len(scene.objects),"| Blender:",bpy.app.version_string)
print("Sources:",[(s["name"],s["xy"],s["strength"]) for s in spillover_sources])
for village in VILLAGES:
    print("Illustrative spillover at village",village["name"],"centre:",
          round(spillover_value(*village["xy"],spillover_sources),4))
    for hi,(hx,hy,_,_) in enumerate(HOUSE_LAYOUT):
        s = spillover_value(village["xy"][0]+hx,village["xy"][1]+hy,spillover_sources)
        inside = [level for level in FIELD_LEVELS if s > level]
        print("   house",village["name"]+str(hi+1),
              "(treated)" if hi+1 in treated_houses[village["name"]] else "         ",
              "S =",format(s,".3f"),
              "| inside S >= "+format(inside[0],"g") if inside else "| outside all levels")
if RENDER_NOW:
    bpy.ops.render.render(scene=scene.name)
