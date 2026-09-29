"""Urban exposure-field counterfactuals -- Blender 5.1.

The original scene construction was render-tested with Blender 4.5.
The lengthscale/variance revision was checked numerically, without rendering.

QUICK START
Scripting > Text > Open this file > Run Script (Alt+P), then F12.
By default this creates six scenes. Choose them with Blender's Scene selector
in the top bar. The baseline is active initially. Existing scenes are preserved.
To build just one: BUILD_ALL_POLICIES=False, then set POLICY below.
All geometry, textures, field values and materials are generated locally;
no add-ons, downloads or external Python packages are required in Blender.

INTERPRETATION
This is one illustrative, smooth, non-negative baseline exposure field a(s)
and five transformations of that SAME field, not six independently drawn fields.
The baseline combines a deterministic mean (two broad hotspots and a background)
with one seeded smooth Gaussian field, represented by finite Fourier features.
FIELD_LENGTHSCALE and FIELD_VARIANCE control its spatial correlation and marginal
variance. The intended covariance is squared-exponential:
  Cov(Z(s),Z(t)) approximately FIELD_VARIANCE * exp(-||s-t||^2/(2*lengthscale^2)).
Conditional on the sampled frequencies, this is a finite Gaussian expansion
with EXACT pointwise variance FIELD_VARIANCE and approximate covariance shape.
Its variance is not the empirical variance across this one map. No per-map
normalisation or clipping is applied. FIELD_VARIANCE=0 disables this component.
Height = FIELD_ZERO_Z + HEIGHT_PER_UNIT * exposure. Colour uses the same fixed
0--EXPOSURE_MAX scale in every scene. The translucent surface is a graph of
exposure over location, not physical terrain, plume height or an outcome effect.
The shallow buildings and streets are geographical context. The field never
intersects them, even at zero. A zero plateau is still a visible surface.

Policies on the natural exposure scale:
  baseline      a(s)
  cap           min(a(s), CAP_LEVEL)
  hinge         max(a(s) - HINGE_DELTA, 0)
  proportionate PROPORTION_REMAINING * a(s)
  local         a(s) - LOCAL_REDUCTION * exp(-||s-LOCAL_CENTRE||^2/(2*LOCAL_SIGMA^2))
  exclusion     0 inside ||s-ZONE_CENTRE|| <= ZONE_RADIUS; a(s) outside

The local policy is an additive smooth translation. Defaults stay positive;
it is NOT clipped at zero, which could introduce another hinge. Invalid sampled
negative values raise an error. Smooth appearance alone does not establish any
statistical support, positivity or Cameron--Martin membership. Exposure values
are illustrative units, not calibrated data or fitted outcomes. This illustrates
counterfactual exposure states; calculating an outcome estimand requires a
response model and the relevant spatial/population averaging measure.

The exclusion policy here means zero EXPOSURE, not a ban on sources followed by
transport/diffusion, and not exclusion of that area from the target population.
The zone boundary is a dense polygon approximating a circle. Separate meshes
represent its one-sided values: there is no smooth interpolation across the jump.
The faint vertical curtain at that boundary is a guide to the discontinuity,
not a continuum of exposure values at the same location.

Rendering follows the attached village diorama: orthographic (18,-18,18) view,
thick soil base, warm sun/cool fill, procedural texture, AgX, studio/alpha output.
Urban geometry, lighting, camera, numerical scale and contour levels are shared.
Labels are optional; numerical scale and pixel anchors help with later annotation.
Changing a policy does not change framing. Surface meshes are finite-resolution
approximations; increase FIELD_ANGULAR_STEPS / reduce FIELD_STEP for close-ups.
No smoothing/subdivision modifier is used that would erase policy plateaux.
"""

import bpy
import bmesh
import json
import math
import os
import random
from functools import lru_cache
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view

# ---------------------------- EDIT THESE ----------------------------
POLICY = "baseline"             # baseline / cap / hinge / proportionate / local / exclusion
BUILD_ALL_POLICIES = True
RENDER_NOW = False               # Render the selected POLICY after building.
RENDER_ALL = False               # Save one PNG per policy; builds all six scenes.
OUTPUT_DIRECTORY = "//exposure_policy_renders"  # Save the .blend first if using //.
OVERWRITE_RENDERS = False        # Otherwise add a numbered suffix to existing files.
WRITE_METADATA = True           # Write matching JSON when RENDER_ALL=True.

SEED = 29
RESOLUTION = (1500, 1150)
SAMPLES = 64                     # 24 for a draft; 128 for a final image.
GROUND_SIZE = (14.8, 12.0)
BACKGROUND = "studio"           # "studio" or "transparent" (keeps board shadow).
BACKDROP_COLOURS = ("C3C8CE", "8A929A")
SUN_STRENGTH = 3.2
SUN_ANGLE = 3.0
SUN_POSITION = (-7, -9, 14)
FILL_RATIO = .33
RIM_STRENGTH = .8
SKY_STRENGTH = .07
COLOUR_LOOK = "AgX - Medium High Contrast"
EXPOSURE = .1                    # Camera exposure, NOT the scientific field.
FRAME_MARGIN = .055

# Baseline a(s) = deterministic mean + seeded smooth Gaussian variation.
# Edit baseline_value(x,y) below to substitute your own analytic field.
FIELD_LENGTHSCALE = 1.2          # Scene units. Smaller -> more, narrower local peaks.
FIELD_VARIANCE = 100.0           # Exposure units squared; 100 means an SD of 10.
FIELD_SEED = 29                  # Fixed realisation in every policy; separate from city SEED.
FIELD_MODES = 64                 # Fourier features; more improves covariance approximation.
BASE_LEVEL = 28.0                # Mean background; raised to keep this realisation positive.
BACKGROUND_WAVE = 3.0
# (centre x, centre y, amplitude, width x, width y); widths are Gaussian SDs.
HOTSPOTS = ((2.0, -.6, 57.0, 2.7, 2.0), (-3.4, 2.5, 34.0, 2.0, 2.35))

CAP_LEVEL = 46.0
HINGE_DELTA = 32.0
PROPORTION_REMAINING = .65       # 35% reduction.
LOCAL_CENTRE = (2.0, -.6)
LOCAL_REDUCTION = 40.0
LOCAL_SIGMA = 1.40
ZONE_CENTRE = (2.0, -.6)
ZONE_RADIUS = 2.10

EXPOSURE_MAX = 120.0             # Common scale with headroom; never re-normalise each panel.
FIELD_ZERO_Z = .58               # Display offset; this elevated plane means ZERO.
HEIGHT_PER_UNIT = .046           # Fixed vertical scale across all policies.
FIELD_INSET = .16
FIELD_STEP = .13                 # Radial mesh spacing; .08 for closer crops.
FIELD_ANGULAR_STEPS = 320
FIELD_OPACITY = .64              # .45 clearer city; .80 clearer surface.
FIELD_COLOURS = ((0.0, "D6EAF1"), (.25, "91C6DD"), (.50, "4597C8"),
                 (.75, "1E69A8"), (1.0, "173E78"))
SHOW_CONTOURS = True
CONTOUR_LEVELS = (20.0, 40.0, 60.0, 80.0, 100.0)
CONTOUR_WIDTH = .012
SHOW_SURFACE_GRID = True         # Sparse curved lines make surface height legible.
SURFACE_GRID_SPACING = 2.0        # Scene-distance units, fixed across all panels.
SHOW_POLICY_BOUNDARY = True      # Cap/hinge kink; exclusion zone edge.
SHOW_ZONE_CURTAIN = True
SHOW_HEIGHT_SCALE = True         # Fixed 0--100 height ruler.
SHOW_TITLE = False               # Prefer adding panel titles in the document.
SHOW_FIELD = True                # False produces the urban base alone.

POLICIES = ("baseline", "cap", "hinge", "proportionate", "local", "exclusion")
TITLES = {"baseline": "Baseline", "cap": "Cap", "hinge": "Hinge",
          "proportionate": "Proportionate reduction", "local": "Smooth local reduction",
          "exclusion": "Exclusion zone"}


@lru_cache(maxsize=8)
def field_modes(seed, count):
    """Freeze frequencies and independent N(0,1) cosine/sine coefficients.

    With M modes, Z(s)=sqrt(v/M)*sum[a_j*cos(w_j.s/ell)+b_j*sin(w_j.s/ell)].
    For fixed frequencies, Var(Z(s))=v since cos^2+sin^2=1. Drawing w_j from
    N(0,I) makes the expected covariance v*exp(-distance^2/(2*ell^2)).
    A finite mode count approximates that covariance, not continuum support.
    """
    rng = random.Random(seed)
    return tuple(tuple(rng.gauss(0,1) for _ in range(4)) for _ in range(count))


@lru_cache(maxsize=131072)
def baseline_value(x, y):
    """One reproducible smooth field, reused at the same locations in all policies.

    Rerunning the script resets this cache. build_scenes() also clears it, so
    calling build_scenes() after changing global settings uses the new field.
    """
    value = BASE_LEVEL + BACKGROUND_WAVE * math.sin(.32*x + .4) * math.cos(.29*y - .2)
    for cx, cy, amplitude, sx, sy in HOTSPOTS:
        value += amplitude * math.exp(-.5 * (((x-cx)/sx)**2 + ((y-cy)/sy)**2))
    if FIELD_VARIANCE:
        variation = 0.0
        for wx,wy,a,b in field_modes(FIELD_SEED,FIELD_MODES):
            phase = (wx*x+wy*y)/FIELD_LENGTHSCALE
            variation += a*math.cos(phase)+b*math.sin(phase)
        value += math.sqrt(FIELD_VARIANCE/FIELD_MODES)*variation
    return value


def policy_value(x, y, policy):
    a = baseline_value(x, y)
    if policy == "baseline":
        return a
    if policy == "cap":
        return min(a, CAP_LEVEL)
    if policy == "hinge":
        return max(a-HINGE_DELTA, 0.0)
    if policy == "proportionate":
        return PROPORTION_REMAINING*a
    if policy == "local":
        d2 = (x-LOCAL_CENTRE[0])**2 + (y-LOCAL_CENTRE[1])**2
        return a-LOCAL_REDUCTION*math.exp(-d2/(2*LOCAL_SIGMA**2))
    if policy == "exclusion":
        d2 = (x-ZONE_CENTRE[0])**2 + (y-ZONE_CENTRE[1])**2
        return 0.0 if d2 <= ZONE_RADIUS**2 else a
    raise ValueError("Unknown POLICY: " + str(policy))


def height(value):
    return FIELD_ZERO_Z + HEIGHT_PER_UNIT*value


def colour(hex_value):
    channels = [int(hex_value[i:i+2], 16)/255 for i in (0, 2, 4)]
    return tuple(c/12.92 if c <= .04045 else ((c+.055)/1.055)**2.4 for c in channels)+(1,)


def collection(name, parent=None):
    coll = bpy.data.collections.new(name)
    if parent:
        parent.children.link(coll)
    return coll


def mesh_object(name, vertices, faces, mat, coll):
    data = bpy.data.meshes.new(name)
    data.from_pydata(vertices, [], faces)
    data.update()
    obj = bpy.data.objects.new(name, data)
    coll.objects.link(obj)
    if mat:
        data.materials.append(mat)
    return obj


CUBE_VERTS = ((-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
              (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1))
CUBE_FACES = ((0,3,2,1), (4,5,6,7), (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7))


def box(name, loc, size, mat, coll, bevel=.018):
    obj = mesh_object(name, [tuple(v[k]*size[k]/2 for k in range(3)) for v in CUBE_VERTS],
                      CUBE_FACES, mat, coll)
    obj.location = loc
    if bevel:
        mod = obj.modifiers.new("Soft edges", "BEVEL")
        mod.width, mod.segments = bevel, 2
    return obj


def strokes(name, lines, mat, coll, width=.012):
    """Many polylines in one curve object, for inexpensive contours."""
    if not lines:
        return None
    data = bpy.data.curves.new(name, "CURVE")
    data.dimensions = "3D"
    data.bevel_depth, data.bevel_resolution = width, 2
    for line in lines:
        if len(line) < 2:
            continue
        spline = data.splines.new("POLY")
        spline.points.add(len(line)-1)
        for point, position in zip(spline.points, line):
            point.co = (*position, 1)
    obj = bpy.data.objects.new(name, data)
    coll.objects.link(obj)
    data.materials.append(mat)
    return obj


def material(name, light, dark=None, scale=5, roughness=.82, bump=.012):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    mat.diffuse_color = colour(light)
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    shader = nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = colour(light)
    shader.inputs["Roughness"].default_value = roughness
    if dark:
        coord = nodes.new("ShaderNodeTexCoord")
        noise = nodes.new("ShaderNodeTexNoise")
        noise.inputs["Scale"].default_value = scale
        noise.inputs["Detail"].default_value = 3
        links.new(coord.outputs["Object"], noise.inputs["Vector"])
        ramp = nodes.new("ShaderNodeValToRGB")
        ramp.color_ramp.elements[0].position = .18
        ramp.color_ramp.elements[0].color = colour(dark)
        ramp.color_ramp.elements[1].position = .82
        ramp.color_ramp.elements[1].color = colour(light)
        links.new(noise.outputs["Fac"], ramp.inputs["Fac"])
        links.new(ramp.outputs["Color"], shader.inputs["Base Color"])
        if bump:
            fine = nodes.new("ShaderNodeTexNoise")
            fine.inputs["Scale"].default_value = 12*scale
            links.new(coord.outputs["Object"], fine.inputs["Vector"])
            relief = nodes.new("ShaderNodeBump")
            relief.inputs["Strength"].default_value = .20
            relief.inputs["Distance"].default_value = bump
            links.new(fine.outputs["Fac"], relief.inputs["Height"])
            links.new(relief.outputs["Normal"], shader.inputs["Normal"])
    return mat


def overlay_material(name, opacity=1.0, constant=None):
    """Fixed value-to-colour mapping; camera-only so field casts no blue shadows."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    mat.diffuse_color = (*colour(constant or "4597C8")[:3], opacity)
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    surface = nodes.new("ShaderNodeBsdfPrincipled")
    surface.inputs["Roughness"].default_value = .43
    surface.inputs["Specular IOR Level"].default_value = .12
    surface.inputs["Emission Strength"].default_value = .85
    if constant:
        surface.inputs["Base Color"].default_value = tuple(c*.18 for c in colour(constant)[:3])+(1,)
        surface.inputs["Emission Color"].default_value = colour(constant)
    else:
        # Derive the colour from the graph's actual world height. This keeps
        # colour and height in exact agreement after triangle interpolation,
        # without depending on named-attribute lookup across Blender versions.
        geometry = nodes.new("ShaderNodeNewGeometry")
        xyz = nodes.new("ShaderNodeSeparateXYZ")
        links.new(geometry.outputs["Position"],xyz.inputs["Vector"])
        offset = nodes.new("ShaderNodeMath")
        offset.operation = "SUBTRACT"
        offset.inputs[1].default_value = FIELD_ZERO_Z
        links.new(xyz.outputs["Z"],offset.inputs[0])
        fraction = nodes.new("ShaderNodeMath")
        fraction.operation = "DIVIDE"
        fraction.inputs[1].default_value = HEIGHT_PER_UNIT*EXPOSURE_MAX
        links.new(offset.outputs[0],fraction.inputs[0])
        ramp = nodes.new("ShaderNodeValToRGB")
        ramp.color_ramp.interpolation = "LINEAR"
        elements = ramp.color_ramp.elements
        for i, (position, hex_value) in enumerate(FIELD_COLOURS):
            element = elements[i] if i < 2 else elements.new(position)
            element.position, element.color = position, colour(hex_value)
        links.new(fraction.outputs[0], ramp.inputs["Fac"])
        dim = nodes.new("ShaderNodeMixRGB")
        dim.blend_type = "MULTIPLY"
        dim.inputs[0].default_value = 1.0
        dim.inputs[2].default_value = (.18,.18,.18,1)
        links.new(ramp.outputs["Color"], dim.inputs[1])
        links.new(dim.outputs["Color"], surface.inputs["Base Color"])
        links.new(ramp.outputs["Color"], surface.inputs["Emission Color"])
    # Modest diffuse shading reveals shape, without refraction or caustics.
    # Colour is a redundant schematic cue; use the ruler for numerical reading.
    rays = nodes.new("ShaderNodeLightPath")
    alpha = nodes.new("ShaderNodeMath")
    alpha.operation = "MULTIPLY"
    alpha.inputs[1].default_value = opacity
    links.new(rays.outputs["Is Camera Ray"], alpha.inputs[0])
    transparent = nodes.new("ShaderNodeBsdfTransparent")
    mix = nodes.new("ShaderNodeMixShader")
    links.new(alpha.outputs[0], mix.inputs[0])
    links.new(transparent.outputs[0], mix.inputs[1])
    links.new(surface.outputs[0], mix.inputs[2])
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(mix.outputs[0], output.inputs["Surface"])
    return mat


def build_city():
    """Shared low-relief city: muted roofs, pavements, road markings and two parks."""
    rng = random.Random(SEED)
    common = collection("Urban base | shared across policies")
    land = collection("01 Ground and soil", common)
    city = collection("02 Shallow urban blocks", common)
    soil = material("Soil | warm earth", "96704C", "634C39", 7, bump=.035)
    subsoil = material("Dark foundation layer", "756C60", "554F47", 10)
    paving = material("Concrete | warm limestone", "D9D4C6", "C2BCAF", 8, bump=.007)
    pavement = material("Pavements | pale stone", "E2DDD0", "CDC7B9", 15, bump=.005)
    road = material("Asphalt | subdued grey", "90989B", "7C858A", 9, bump=.006)
    paint = material("Road markings | ivory", "E4E0D1")
    walls = material("Buildings | warm plaster", "E3D8C1", "C9C0AE", 8, bump=.004)
    roofs = [material("Roof | " + h, h, l, 10, bump=.006) for h,l in
             (("BDBCAF","A8AA9F"), ("C3B7A4","AEA18E"), ("B49C86","9E8976"))]
    glass = material("Roof lights and windows", "74898B", roughness=.5)
    green = material("Pocket parks | muted sage", "A7B18D", "8D9D77", 6, bump=.006)
    leaf = material("Tree crowns | muted olive", "9EAF86", "7E956C", 8, bump=.007)
    gx, gy = GROUND_SIZE
    box("Exposed earth slab", (0,0,-.43), (gx,gy,.70), soil, land, .15)
    box("Foundation seam", (0,0,-.10), (gx-.025,gy-.025,.15), subsoil, land, .045)
    box("City ground", (0,0,-.012), (gx-.055,gy-.055,.05), paving, land, .035)
    xcuts, ycuts = (-7.0,-3.5,0,3.5,7.0), (-5.6,-1.85,1.85,5.6)
    for x in xcuts[1:-1]:
        box("North south street", (x,0,.021), (.57,gy-.25,.016), road, city, .006)
        for k in range(-11,12):
            y = k*.48
            if all(abs(y-yc) > .46 for yc in ycuts[1:-1]):
                box("Lane dash", (x,y,.032), (.022,.16,.004), paint, city, 0)
    for y in ycuts[1:-1]:
        box("East west street", (0,y,.023), (gx-.25,.57,.016), road, city, .006)
        for k in range(-14,15):
            x = k*.48
            if all(abs(x-xc) > .48 for xc in xcuts[1:-1]):
                box("Lane dash", (x,y,.034), (.16,.022,.004), paint, city, 0)
    for ix in range(4):
        for iy in range(3):
            x0,x1,y0,y1 = xcuts[ix],xcuts[ix+1],ycuts[iy],ycuts[iy+1]
            cx,cy = (x0+x1)/2,(y0+y1)/2
            sx,sy = x1-x0-.83,y1-y0-.83
            box("Raised pavement", (cx,cy,.053), (sx+.13,sy+.13,.065), pavement, city, .035)
            if (ix,iy) in ((0,0),(3,2)):
                box("Pocket park", (cx,cy,.093), (sx-.16,sy-.16,.022), green, city, .025)
                box("Park path", (cx,cy,.109), (.19,sy-.20,.012), paving, city, .005)
                for tx,ty in ((-.68,-.7),(.63,-.62),(-.62,.62),(.69,.62)):
                    box("Tree trunk", (cx+tx,cy+ty,.16), (.045,.045,.13), subsoil, city, .008)
                    data = bpy.data.meshes.new("Low tree crown")
                    bm = bmesh.new()
                    bmesh.ops.create_icosphere(bm, subdivisions=2, radius=.18)
                    bm.to_mesh(data)
                    bm.free()
                    obj = bpy.data.objects.new("Low tree crown", data)
                    city.objects.link(obj)
                    obj.location = (cx+tx,cy+ty,.30)
                    obj.scale = (1,1,.75)
                    data.materials.append(leaf)
                continue
            for bx in (-1,1):
                for by in (-1,1):
                    x,y = cx+bx*sx*.245,cy+by*sy*.245
                    wx,wy = sx*.40,sy*.39
                    h = rng.uniform(.15,.26)
                    top = .087+h
                    roof = roofs[rng.randrange(len(roofs))]
                    box("Urban building", (x,y,.087+h/2), (wx,wy,h), walls, city, .022)
                    box("Flat roof parapet", (x,y,top+.014), (wx+.045,wy+.045,.040), walls, city, .009)
                    box("Recessed roof", (x,y,top+.037), (wx-.07,wy-.07,.012), roof, city, .008)
                    box("Roof light", (x+.16,y,top+.049), (wx*.25,wy*.40,.014), glass, city, .006)
                    for dx in (-.25,.25):
                        box("Facade window", (x+dx,y-wy/2-.003,top-.065), (.13,.006,.055), glass, city, .001)
    return common


def planar_grid():
    """One grid for all panels, with an explicit ring at the exclusion boundary.

    Ray endpoints include the rectangle's four corners, so the footprint is
    exactly rectangular. The circle uses straight segments between ring points.
    Returns xy vertices, inside faces, outside faces, ring indices, perimeter.
    """
    hx,hy = GROUND_SIZE[0]/2-FIELD_INSET, GROUND_SIZE[1]/2-FIELD_INSET
    cx,cy = ZONE_CENTRE
    if not (abs(cx)+ZONE_RADIUS < hx and abs(cy)+ZONE_RADIUS < hy):
        raise ValueError("The entire exclusion zone must lie inside the board.")
    angles = [2*math.pi*i/FIELD_ANGULAR_STEPS for i in range(FIELD_ANGULAR_STEPS)]
    angles += [math.atan2(y-cy,x-cx) % (2*math.pi) for x in (-hx,hx) for y in (-hy,hy)]
    angles = sorted(set(round(a,12) for a in angles))
    n = len(angles)
    extents = []
    for a in angles:
        dx,dy = math.cos(a), math.sin(a)
        tx = ((hx if dx > 0 else -hx)-cx)/dx if abs(dx)>1e-10 else 1e10
        ty = ((hy if dy > 0 else -hy)-cy)/dy if abs(dy)>1e-10 else 1e10
        extents.append(min(tx,ty))
    ni = max(4, math.ceil(ZONE_RADIUS/FIELD_STEP))
    no = max(4, math.ceil((max(extents)-ZONE_RADIUS)/FIELD_STEP))
    xy = [(cx,cy)]
    rings = []
    for j in range(1,ni+no+1):
        ring = []
        for i,a in enumerate(angles):
            r = ZONE_RADIUS*j/ni if j <= ni else ZONE_RADIUS+(extents[i]-ZONE_RADIUS)*(j-ni)/no
            ring.append(len(xy))
            xy.append((cx+r*math.cos(a), cy+r*math.sin(a)))
        rings.append(ring)
    inside = [(0,rings[0][i],rings[0][(i+1)%n]) for i in range(n)]
    outside = []
    for j in range(1,len(rings)):
        target = inside if j < ni else outside
        for i in range(n):
            a,b,c,d = rings[j-1][i],rings[j-1][(i+1)%n],rings[j][(i+1)%n],rings[j][i]
            target.extend(((a,d,c),(a,c,b)))
    return xy,inside,outside,rings[ni-1],rings[-1]


def contour_segments(xy, faces, values, level, z_value):
    lines = []
    for face in faces:
        hits = []
        for a,b in zip(face,face[1:]+face[:1]):
            va,vb = values[a],values[b]
            if (va < level <= vb) or (vb < level <= va):
                t = (level-va)/(vb-va)
                hits.append((xy[a][0]+t*(xy[b][0]-xy[a][0]),
                             xy[a][1]+t*(xy[b][1]-xy[a][1]), z_value))
        if len(hits) == 2 and (Vector(hits[0])-Vector(hits[1])).length > 1e-7:
            lines.append(hits)
    return lines


def field_mesh(name, xy, faces, values, mat, coll):
    # Compact unused vertices, important for the two sides of the exclusion jump.
    used = sorted({i for f in faces for i in f})
    remap = {old:new for new,old in enumerate(used)}
    obj = mesh_object(name, [(xy[i][0],xy[i][1],height(values[i])) for i in used],
                      [tuple(remap[i] for i in f) for f in faces], mat, coll)
    attribute = obj.data.attributes.new("exposure_fraction", "FLOAT", "POINT")
    attribute.data.foreach_set("value", [values[i]/EXPOSURE_MAX for i in used])
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    return obj


def surface_grid_lines(policy):
    """Curves along fixed x/y transects; split at the exclusion jump exactly."""
    hx,hy = GROUND_SIZE[0]/2-FIELD_INSET, GROUND_SIZE[1]/2-FIELD_INSET
    lines = []
    for direction in (0,1):
        fixed_extent,travel_extent = (hx,hy) if direction == 0 else (hy,hx)
        centre_fixed,centre_travel = ZONE_CENTRE if direction == 0 else ZONE_CENTRE[::-1]
        count = math.floor(fixed_extent/SURFACE_GRID_SPACING)
        for k in range(-count,count+1):
            fixed = k*SURFACE_GRID_SPACING
            cuts = [-travel_extent,travel_extent]
            remaining = ZONE_RADIUS**2-(fixed-centre_fixed)**2
            if policy == "exclusion" and remaining > 0:
                cuts += [centre_travel-math.sqrt(remaining),centre_travel+math.sqrt(remaining)]
            cuts.sort()
            for lo,hi in zip(cuts,cuts[1:]):
                mid = (lo+hi)/2
                inner = policy == "exclusion" and (fixed-centre_fixed)**2+(mid-centre_travel)**2 < ZONE_RADIUS**2
                steps = max(2,math.ceil((hi-lo)/.055))
                line = []
                for i in range(steps+1):
                    travel = lo+(hi-lo)*i/steps
                    x,y = (fixed,travel) if direction == 0 else (travel,fixed)
                    value = (0.0 if inner else baseline_value(x,y)) if policy == "exclusion" else policy_value(x,y,policy)
                    line.append((x,y,height(value)+.017))
                lines.append(line)
    return lines


def build_field(policy, grid, parent, field_mat, line_mat, edge_mat):
    coll = collection("03 Exposure | " + policy, parent)
    xy,inside,outside,ring,perimeter = grid
    faces = inside+outside
    original = [baseline_value(*p) for p in xy]
    values = [policy_value(*p,policy) for p in xy]
    # Explicit one-sided limits: do not let floating-point circle tests place
    # boundary vertices on the wrong side, or interpolate across the jump.
    if policy == "exclusion":
        inner_values = [0.0]*len(xy)
        field_mesh("Zero exposure inside zone", xy, inside, inner_values, field_mat, coll)
        field_mesh("Unchanged field outside zone", xy, outside, original, field_mat, coll)
        contour_faces, contour_values = outside, original
    else:
        field_mesh("Exposure surface | " + policy, xy, faces, values, field_mat, coll)
        contour_faces, contour_values = faces, values
    if SHOW_SURFACE_GRID:
        grid_mat = overlay_material("Surface grid | " + policy,.52,"2C6082")
        strokes("Surface coordinate grid",surface_grid_lines(policy),grid_mat,coll,.010)
    if SHOW_CONTOURS:
        lines = []
        for level in CONTOUR_LEVELS:
            lines += contour_segments(xy,contour_faces,contour_values,level,height(level)+.009)
        strokes("Fixed exposure contours",lines,line_mat,coll,CONTOUR_WIDTH)
    if SHOW_POLICY_BOUNDARY and policy in ("cap","hinge"):
        level = CAP_LEVEL if policy == "cap" else HINGE_DELTA
        z = height(CAP_LEVEL if policy == "cap" else 0)+.012
        lines = contour_segments(xy,faces,original,level,z)
        strokes("Policy plateau boundary",lines,edge_mat,coll,CONTOUR_WIDTH*1.5)
    # Show a thin rim, not a solid vertical slab, along the field's outer extent.
    rim = [(xy[i][0],xy[i][1],height(values[i])+.005) for i in perimeter]
    strokes("Field perimeter", [rim+[rim[0]]], line_mat, coll, .011)
    if policy == "exclusion":
        low = [(xy[i][0],xy[i][1],height(0)+.010) for i in ring]
        high = [(xy[i][0],xy[i][1],height(original[i])+.010) for i in ring]
        if SHOW_POLICY_BOUNDARY:
            strokes("Exclusion boundary | both one-sided values",[low+[low[0]],high+[high[0]]],edge_mat,coll,.020)
            ground = [(x,y,.12) for x,y,z in low]
            strokes("Exclusion footprint on city",[ground+[ground[0]]],edge_mat,coll,.017)
        if SHOW_ZONE_CURTAIN:
            n = len(ring)
            mesh_object("Discontinuity guide | not an interpolated field",low+high,
                        [(i,(i+1)%n,(i+1)%n+n,i+n) for i in range(n)],
                        overlay_material("Exclusion boundary curtain",.13,"4597C8"),coll)
    return values


def aim(obj, target=(0,0,0)):
    obj.rotation_euler = (Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()


def text_object(body, position, size, mat, coll, camera, align="LEFT"):
    data = bpy.data.curves.new(body, "FONT")
    data.body, data.size, data.align_x = body, size, align
    obj = bpy.data.objects.new(body, data)
    coll.objects.link(obj)
    obj.location = position
    obj.rotation_euler = camera.rotation_euler
    data.materials.append(mat)
    return obj


def build_studio(grid):
    studio = collection("04 Studio and common scale")
    floor_mat = material("Backdrop | cool grey", BACKDROP_COLOURS[0])
    nodes,links = floor_mat.node_tree.nodes,floor_mat.node_tree.links
    geometry = nodes.new("ShaderNodeNewGeometry")
    planar = nodes.new("ShaderNodeVectorMath")
    planar.operation = "MULTIPLY"
    planar.inputs[1].default_value = (1,1,0)
    links.new(geometry.outputs["Position"],planar.inputs[0])
    distance = nodes.new("ShaderNodeVectorMath")
    distance.operation = "LENGTH"
    links.new(planar.outputs[0],distance.inputs[0])
    mapping = nodes.new("ShaderNodeMapRange")
    reach = math.hypot(*GROUND_SIZE)/2
    mapping.inputs["From Min"].default_value,mapping.inputs["From Max"].default_value = reach*.75,reach*2.1
    links.new(distance.outputs["Value"],mapping.inputs["Value"])
    ramp = nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.interpolation = "EASE"
    for elem,h in zip(ramp.color_ramp.elements,BACKDROP_COLOURS):
        elem.color = colour(h)
    links.new(mapping.outputs["Result"],ramp.inputs["Fac"])
    links.new(ramp.outputs["Color"],nodes.get("Principled BSDF").inputs["Base Color"])
    floor = box("Studio floor",(0,0,-.83),(200,200,.10),floor_mat,studio,0)
    floor.is_shadow_catcher = BACKGROUND == "transparent"
    sun_vector = Vector(SUN_POSITION).normalized()
    data = bpy.data.lights.new("Warm key sun","SUN")
    data.energy,data.angle,data.color = SUN_STRENGTH,math.radians(SUN_ANGLE),colour("FFE8CC")[:3]
    sun = bpy.data.objects.new("Warm key sun",data)
    studio.objects.link(sun)
    sun.location = sun_vector*15
    aim(sun)
    for name,loc,energy,size,tint in (("Cool fill",(7,1,9),SUN_STRENGTH*FILL_RATIO*260,9,"DCE8FF"),
                                     ("Cool rim",(-10,11,7),RIM_STRENGTH*650,5,"D2E2FF")):
        data = bpy.data.lights.new(name,"AREA")
        data.energy,data.shape,data.size,data.color = energy,"DISK",size,colour(tint)[:3]
        obj = bpy.data.objects.new(name,data)
        studio.objects.link(obj)
        obj.location = loc
        aim(obj)
    world = bpy.data.worlds.new("Diorama daylight | built-in sky")
    world.use_nodes = True
    sky = world.node_tree.nodes.new("ShaderNodeTexSky")
    available = [item.identifier for item in sky.bl_rna.properties["sky_type"].enum_items]
    for kind in ("MULTIPLE_SCATTERING","NISHITA","SINGLE_SCATTERING","HOSEK_WILKIE"):
        if kind in available:
            sky.sky_type = kind
            break
    if hasattr(sky,"sun_disc"):
        sky.sun_disc = False
    if hasattr(sky,"sun_elevation"):
        sky.sun_elevation = math.asin(sun_vector.z)
        sky.sun_rotation = math.atan2(-sun_vector.x,sun_vector.y) % (2*math.pi)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = SKY_STRENGTH
    world.node_tree.links.new(sky.outputs["Color"],world.node_tree.nodes["Background"].inputs["Color"])
    data = bpy.data.cameras.new("Isometric camera | identical in every panel")
    data.type,data.lens,data.clip_end = "ORTHO",50,500
    camera = bpy.data.objects.new("Isometric camera",data)
    studio.objects.link(camera)
    target = Vector((0,0,.15))
    camera.location = target+Vector((18,-18,18))
    aim(camera,target)
    rotation = camera.rotation_euler.to_matrix()
    right,up = rotation@Vector((1,0,0)),rotation@Vector((0,1,0))
    # Fit the board plus the BASELINE surface and full ruler, never the current
    # policy. All policies lie between baseline and zero, so they fit this same
    # frame without reserving an empty box above every corner of the board.
    hx,hy = GROUND_SIZE[0]/2,GROUND_SIZE[1]/2
    points = [Vector((x,y,z)) for x in (-hx,hx) for y in (-hy,hy) for z in (-.80,.06)]
    points += [Vector((x,y,height(baseline_value(x,y)))) for x,y in grid[0]]
    for value in (0,EXPOSURE_MAX):
        p = Vector((hx-.18,hy-.28,height(value)))
        points.extend((p+right*.85,p+up*.60-right*1.70))
    us,vs = [p.dot(right) for p in points],[p.dot(up) for p in points]
    width,vertical = max(us)-min(us),max(vs)-min(vs)
    x,y = RESOLUTION
    fit = max(width,vertical*x/y) if x>=y else max(width*y/x,vertical)
    data.ortho_scale = fit/(1-2*FRAME_MARGIN)
    camera.location += right*((max(us)+min(us))/2-camera.location.dot(right))
    camera.location += up*((max(vs)+min(vs))/2-camera.location.dot(up))
    if SHOW_HEIGHT_SCALE:
        ink = overlay_material("Scale and annotations",1,"374A56")
        x,y = hx-.18,hy-.28
        axis = [(x,y,height(0)),(x,y,height(EXPOSURE_MAX))]
        lines = [axis]
        for value in (0,.25*EXPOSURE_MAX,.5*EXPOSURE_MAX,.75*EXPOSURE_MAX,EXPOSURE_MAX):
            p = Vector((x,y,height(value)))
            lines.append([tuple(p-right*.08),tuple(p+right*.08)])
            text_object(format(value,"g"),p+right*.16-up*.07,.24,ink,studio,camera)
        strokes("Exposure height ruler",lines,ink,studio,.012)
        p = Vector((x,y,height(EXPOSURE_MAX)))+up*.24-right*.30
        text_object("Exposure",p,.25,ink,studio,camera,align="RIGHT")
    return studio,camera,world


def validate(grid):
    if POLICY not in POLICIES:
        raise ValueError("POLICY must be one of " + str(POLICIES))
    if not (0 <= PROPORTION_REMAINING <= 1 and CAP_LEVEL >= 0 and HINGE_DELTA >= 0):
        raise ValueError("Invalid reduction parameters.")
    if LOCAL_SIGMA <= 0 or ZONE_RADIUS <= 0 or LOCAL_REDUCTION < 0:
        raise ValueError("Widths must be positive and local reduction non-negative.")
    if EXPOSURE_MAX <= 0 or HEIGHT_PER_UNIT <= 0 or not (0 <= FIELD_OPACITY <= 1):
        raise ValueError("Invalid field display settings.")
    if SURFACE_GRID_SPACING <= 0:
        raise ValueError("SURFACE_GRID_SPACING must be positive.")
    if FIELD_ZERO_Z <= .49:
        raise ValueError("FIELD_ZERO_Z must be above the shallow buildings and trees (>.49).")
    values = [baseline_value(*p) for p in grid[0]]
    if min(values) < 0:
        raise ValueError(f"Sampled baseline minimum {min(values):.2f} is negative. Raise BASE_LEVEL or reduce FIELD_VARIANCE; no clipping was applied.")
    if max(values) > EXPOSURE_MAX:
        raise ValueError(f"Sampled baseline maximum {max(values):.2f} exceeds EXPOSURE_MAX. Increase that common scale for all panels.")
    reduced = [policy_value(*p,"local") for p in grid[0]]
    if min(reduced) < 0:
        raise ValueError("Local reduction creates negative exposure. Reduce LOCAL_REDUCTION; no clipping was applied.")
    if BACKGROUND not in ("studio","transparent"):
        raise ValueError("BACKGROUND must be studio or transparent.")


def policy_formula(policy):
    return {"baseline":"a(s)", "cap":f"min(a(s), {CAP_LEVEL:g})",
            "hinge":f"max(a(s) - {HINGE_DELTA:g}, 0)",
            "proportionate":f"{PROPORTION_REMAINING:g} * a(s)",
            "local":f"a(s) - {LOCAL_REDUCTION:g} * exp(-distance(s, {LOCAL_CENTRE})^2 / (2 * {LOCAL_SIGMA:g}^2))",
            "exclusion":f"0 when distance(s, {ZONE_CENTRE}) <= {ZONE_RADIUS:g}; a(s) otherwise"}[policy]


def scene_metadata(scene,policy,values,camera):
    scene.view_layers[0].update()
    def pixel(p):
        v = world_to_camera_view(scene,camera,Vector(p))
        return [round(v.x*RESOLUTION[0],1),round((1-v.y)*RESOLUTION[1],1)]
    points = {"local centre at zero":pixel((*LOCAL_CENTRE,height(0))),
              "zone centre at zero":pixel((*ZONE_CENTRE,height(0))),
              "origin at zero":pixel((0,0,height(0)))}
    record = {"policy":policy,"formula":policy_formula(policy),"units":"illustrative exposure units",
              "baseline":"shared deterministic mean plus one seeded finite-Fourier Gaussian field; approximate squared-exponential covariance",
              "baseline_settings":{"base":BASE_LEVEL,"wave":BACKGROUND_WAVE,"hotspots":HOTSPOTS,
                                   "lengthscale":FIELD_LENGTHSCALE,"marginal_variance":FIELD_VARIANCE,
                                   "field_seed":FIELD_SEED,"fourier_modes":FIELD_MODES,
                                   "covariance_target":"variance * exp(-distance^2 / (2 * lengthscale^2))"},
              "parameters":{"cap":CAP_LEVEL,"hinge_delta":HINGE_DELTA,"proportion_remaining":PROPORTION_REMAINING,
                            "local_centre":LOCAL_CENTRE,"local_reduction":LOCAL_REDUCTION,"local_sigma":LOCAL_SIGMA,
                            "zone_centre":ZONE_CENTRE,"zone_radius":ZONE_RADIUS},
              "common_scale":[0,EXPOSURE_MAX],"field_zero_z":FIELD_ZERO_Z,"height_per_unit":HEIGHT_PER_UNIT,
              "colour_stops":FIELD_COLOURS,"opacity":FIELD_OPACITY,"contour_levels":CONTOUR_LEVELS,
              "sampled_range":[min(values),max(values)],"image_size":list(RESOLUTION),
              "pixel_origin":"top-left; y down","label_anchors":points,
              "camera_location":list(camera.location),"camera_rotation":list(camera.rotation_euler),
              "camera_ortho_scale":camera.data.ortho_scale,"blender":bpy.app.version_string}
    scene["policy"] = policy
    scene["field_metadata"] = json.dumps(record)
    return record


def available_output(directory,stem):
    path = os.path.join(directory,stem)
    if OVERWRITE_RENDERS:
        return path
    i = 1
    while os.path.exists(path+".png") or os.path.exists(path+".json"):
        path = os.path.join(directory,stem+f"_{i:03d}")
        i += 1
    return path


def build_scenes():
    """Public entry point; creates scenes and returns them indexed by policy."""
    if not (math.isfinite(FIELD_LENGTHSCALE) and FIELD_LENGTHSCALE > 0):
        raise ValueError("FIELD_LENGTHSCALE must be finite and positive.")
    if not (math.isfinite(FIELD_VARIANCE) and FIELD_VARIANCE >= 0):
        raise ValueError("FIELD_VARIANCE must be finite and non-negative.")
    if not isinstance(FIELD_MODES,int) or FIELD_MODES < 1:
        raise ValueError("FIELD_MODES must be a positive integer.")
    baseline_value.cache_clear()
    if FIELD_STEP <= 0 or FIELD_ANGULAR_STEPS < 32 or ZONE_RADIUS <= 0:
        raise ValueError("Use FIELD_STEP>0, FIELD_ANGULAR_STEPS>=32 and ZONE_RADIUS>0.")
    if LOCAL_SIGMA <= 0 or any(sx <= 0 or sy <= 0 for _,_,_,sx,sy in HOTSPOTS):
        raise ValueError("Gaussian widths must be positive.")
    grid = planar_grid()
    validate(grid)
    city = build_city()
    studio,camera,world = build_studio(grid)
    field_mat = overlay_material("Exposure | fixed common scale",FIELD_OPACITY)
    line_mat = overlay_material("Exposure contours",.75,"245B82")
    edge_mat = overlay_material("Policy boundary",.95,"204D6D")
    scenes = {}
    chosen = POLICIES if BUILD_ALL_POLICIES or RENDER_ALL else (POLICY,)
    for number,policy in enumerate(POLICIES):
        if policy not in chosen:
            continue
        scene = bpy.data.scenes.new(f"Field {number} | {TITLES[policy]}")
        scenes[policy] = scene
        scene.collection.children.link(city)
        scene.collection.children.link(studio)
        scene.camera,scene.world = camera,world
        scene.render.engine = "CYCLES"
        scene.cycles.device = "CPU"
        scene.cycles.samples,scene.cycles.seed = SAMPLES,SEED
        scene.cycles.use_denoising = True
        scene.cycles.transparent_max_bounces = 24
        scene.render.resolution_x,scene.render.resolution_y = RESOLUTION
        scene.render.resolution_percentage = 100
        scene.render.image_settings.file_format = "PNG"
        scene.render.image_settings.color_mode = "RGBA"
        scene.render.film_transparent = BACKGROUND == "transparent"
        scene.view_settings.view_transform = "AgX"
        try:
            scene.view_settings.look = COLOUR_LOOK
        except TypeError:
            print("Using default AgX look; requested look unavailable:",COLOUR_LOOK)
        scene.view_settings.exposure = EXPOSURE
        if SHOW_FIELD:
            values = build_field(policy,grid,scene.collection,field_mat,line_mat,edge_mat)
        else:
            values = [policy_value(*p,policy) for p in grid[0]]
        if SHOW_TITLE:
            labels = collection("05 Optional title",scene.collection)
            rotation = camera.rotation_euler.to_matrix()
            up = rotation@Vector((0,1,0))
            visible_height = camera.data.ortho_scale*min(1,RESOLUTION[1]/RESOLUTION[0])
            forward = rotation@Vector((0,0,-1))
            location = camera.location+forward*25+up*(visible_height/2-.55)
            text_object(TITLES[policy],location,.35,edge_mat,labels,camera,"CENTER")
        record = scene_metadata(scene,policy,values,camera)
        print(scene.name,"|",record["formula"],"| sampled range",[round(v,3) for v in record["sampled_range"]])
    if bpy.context.window:
        bpy.context.window.scene = scenes[POLICY]
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "VIEW_3D":
                space = area.spaces.active
                space.shading.type = "MATERIAL"
                space.shading.use_scene_world = True
                space.shading.use_scene_lights = True
                space.overlay.show_overlays = False
                space.region_3d.view_perspective = "CAMERA"
                space.region_3d.view_camera_zoom = 0
    return scenes


def main():
    scenes = build_scenes()
    if RENDER_ALL:
        if OUTPUT_DIRECTORY.startswith("//") and not bpy.data.filepath:
            raise ValueError("Scenes built. Save the .blend first, or set OUTPUT_DIRECTORY to an absolute folder.")
        directory = bpy.path.abspath(OUTPUT_DIRECTORY)
        os.makedirs(directory,exist_ok=True)
        for i,policy in enumerate(POLICIES):
            scene = scenes[policy]
            stem = available_output(directory,f"{i:02d}_{policy}")
            scene.render.filepath = stem+".png"
            bpy.ops.render.render(scene=scene.name,write_still=True)
            if WRITE_METADATA:
                with open(stem+".json","w",encoding="utf-8") as handle:
                    json.dump(json.loads(scene["field_metadata"]),handle,indent=2)
        if bpy.context.window:
            bpy.context.window.scene = scenes[POLICY]
    elif RENDER_NOW:
        bpy.ops.render.render(scene=scenes[POLICY].name)
    print("Urban exposure scenes ready. Choose a scene in the top bar, then press F12.")
    print("Height and colours share a fixed scale. All scenes use the same city, lights and camera.")
    print("Save the .blend with File > Save As. Existing scenes have been preserved.")
    return scenes


if __name__ == "__main__":
    main()
