import argparse
import json
import vtk
import numpy
from numpy import sqrt, sin, cos, pi, arccos, arcsin
import numpy as np
import vtk
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkUnsignedCharArray
import cppyy
from cppyy.gbl.std import vector
cppyy.include("quaternion.h")
quaternion=cppyy.gbl.SPPARKS_NS.quaternion

import bor_variants as bor
# Allowable variants ['V1', 'V2', ..., 'V12']
# VALID_VARIANTS = [f'V{i}' for i in range(1, 13)]
VALID_VARIANTS=bor.VALID_VARIANTS

def get_tube_filter(polydata):
    # This is work around on Mac M3 Max that won't render wireframe.
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(polydata)
    # Adjust tube thickness as needed
    tubeFilter.SetRadius(0.015)
    # Number of sides for the tube (higher = smoother)
    # In this case more is not needed.
    tubeFilter.SetNumberOfSides(1)
    tubeFilter.Update()
    return tubeFilter

def get_plane_actor(center, normal, size=1.0):
    """
    Creates a VTK actor for a plane centered at a given point with a specified
    normal.

    Args:
        center (np.array): A 3-element numpy array representing the
        center of the plane.
        normal (np.array): A 3-element numpy array representing the
        normal vector of the plane.
        size (float): The approximate side length of the square plane.
        The actual size will be adjusted by VTK's plane source.

    Returns:
        vtkActor: A VTK actor representing the plane.
    """
    # Normalize the normal vector
    normal = normal / np.linalg.norm(normal)

    # Create a plane source
    planeSource = vtk.vtkPlaneSource()
    planeSource.SetOrigin(-size / 2, -size / 2, 0)
    planeSource.SetPoint1(size / 2, -size / 2, 0)
    planeSource.SetPoint2(-size / 2, size / 2, 0)

    # Compute a transformation matrix to orient and position the plane
    # The default plane source is in the XY plane. We need to rotate it
    # so its normal aligns with the given normal and then translate it.

    # 1. Determine the rotation to align the z-axis with the normal
    # The normal of the plane source is initially (0, 0, 1).
    # We want to rotate it so (0, 0, 1) goes to normal.
    z_axis = np.array([0, 0, 1])
    rotation_axis = numpy.cross(z_axis, normal)
    rotation_angle = numpy.degrees(np.arccos(np.dot(z_axis, normal)))

    transform = vtk.vtkTransform()

    # Avoid division by zero if vectors are parallel 
    if np.linalg.norm(rotation_axis) > 1e-6: 
        transform.RotateWXYZ(rotation_angle, rotation_axis[0],
                             rotation_axis[1], rotation_axis[2])

    # 2. Translate the plane so its center is at the desired center
    # The plane source is initially centered at (0, 0, 0).
    transform.Translate(center[0], center[1], center[2])

    # Apply the transformation to the plane source's output
    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetInputConnection(planeSource.GetOutputPort())
    transformFilter.SetTransform(transform)
    transformFilter.Update()

    # Create a mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transformFilter.GetOutputPort())

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # Red color for the plane
    actor.GetProperty().SetColor(1.0, 0.0, 0.0)  
    # Semi-transparent
    actor.GetProperty().SetOpacity(1.0)         
    # Ensure good lighting interaction
    actor.GetProperty().SetDiffuseColor(1.0, 0.0, 0.0) 
    actor.GetProperty().SetAmbientColor(0.2, 0.0, 0.0)
    actor.GetProperty().SetSpecular(0.5)
    actor.GetProperty().SetSpecularPower(50)
    return actor

def get_hexagonal_prism_actor(orient=None,center=(0,0,0),side=1.0,height=2.0):
    # orient: quaternion, numpy array; used to orient prism
    q=orient
    # print("%8.6f, %8.6f, %8.6f, %8.6f"%(q[0], q[1], q[2], q[3]))
    s = side
    h = height
 
    # prism centered on this point
    x0=center[0]
    y0=center[1]
    z0=center[2]

    points = vtk.vtkPoints()

    # Bottom hexagon vertices (z=-h/2)
    z=-h/2.
    for i in range(6):
        angle = 2 * pi * i / 6
        x = s * cos(angle)
        y = s * sin(angle)
        xyz=vector['double']((x,y,z))
        p=cppyy.gbl.SPPARKS_NS.quaternion.rotate_vector(q,xyz);
        points.InsertNextPoint(x0+p[0],y0+p[1],z0+p[2])

    # Top hexagon vertices (z=h/2)
    z=h/2.
    for i in range(6):
        angle = 2 * pi * i / 6
        x = s * cos(angle)
        y = s * sin(angle)
        xyz=vector['double']((x,y,z))
        p=cppyy.gbl.SPPARKS_NS.quaternion.rotate_vector(q,xyz);
        points.InsertNextPoint(x0+p[0],y0+p[1],z0+p[2])

    # Define the lines (edges) of the hexagonal prism
    lines = vtk.vtkCellArray()
    colors = vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    namedColors=vtkNamedColors()

    # Bottom hexagon edges
    for i in range(6):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, i)
        line.GetPointIds().SetId(1, (i + 1) % 6)
        lines.InsertNextCell(line)
        colors.InsertNextTypedTuple(namedColors.GetColor3ub("Black"))

    # Top hexagon edges (indices 6 to 11 for top vertices)
    for i in range(6):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 6 + i)
        line.GetPointIds().SetId(1, 6 + ((i + 1) % 6))
        lines.InsertNextCell(line)
        colors.InsertNextTypedTuple(namedColors.GetColor3ub("Black"))

    # Hack in extra line along [11-20]
    # line = vtk.vtkLine()
    # line.GetPointIds().SetId(0, 6)
    # line.GetPointIds().SetId(1, 9)
    # lines.InsertNextCell(line)
    # colors.InsertNextTypedTuple(namedColors.GetColor3ub("Red"))
    # END Hack in extra line along a1

    # Vertical edges connecting bottom and top hexagons
    for i in range(6):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, i)
        line.GetPointIds().SetId(1, 6 + i)
        lines.InsertNextCell(line)
        colors.InsertNextTypedTuple(namedColors.GetColor3ub("Black"))

    # # Create a sphere source for the glyph
    ## JAK: Save -- show how to create a glyph actor
    # sphere_source = vtk.vtkSphereSource()
    # glyph_radius=.05*h
    # sphere_source.SetRadius(glyph_radius)
    # sphere_source.SetThetaResolution(16)
    # sphere_source.SetPhiResolution(16)
    # sphere_source.Update() # Important to update the source
    #
    # # Create vtkGlyph3D for the sphere glyph
    # # We want the glyph only at the end_point.
    # # So, we create a vtkPoints object containing only the end_point.
    # # glyph point is 1st point
    # glyph_point_id=6
    # glyph_points = vtk.vtkPoints()
    # glyph_points.InsertNextPoint(points.GetPoint(glyph_point_id))
    #
    # glyph_poly_data = vtk.vtkPolyData()
    # glyph_poly_data.SetPoints(glyph_points)
    #
    # glyph3d = vtk.vtkGlyph3D()
    # glyph3d.SetSourceConnection(sphere_source.GetOutputPort())
    # glyph3d.SetInputData(glyph_poly_data) # Use the poly data with just the end point
    # glyph3d.Update()
    #
    # # Glyph mapper and actor
    # glyph_mapper = vtk.vtkPolyDataMapper()
    # glyph_mapper.SetInputConnection(glyph3d.GetOutputPort())
    #
    # glyph_actor = vtk.vtkActor()
    # glyph_actor.SetMapper(glyph_mapper)
    # glyph_actor.GetProperty().SetColor(namedColors.GetColor3ub("Blue"))
    #

    # Create a PolyData object to store the points and lines
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetCellData().SetScalars(colors)

    # Create a prism mapper
    mapper = vtk.vtkPolyDataMapper()
    # mapper.SetInputData(polydata)
    tube_filter=get_tube_filter(polydata)
    mapper.SetInputConnection(tube_filter.GetOutputPort())

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetLineWidth(4.0) # Make the lines a bit thicker
    return actor

def get_cube_actor(orient=vector['double']((1,0,0,0)),center=(0,0,0),edge=2.0):

    # orient: quaternion, numpy array; used to orient cube
    q=orient
    # length of cube edge
    s = edge

    # prism centered on this point
    x0=center[0]
    y0=center[1]
    z0=center[2]

    local_points = [
        (-s/2, -s/2, -s/2),  # 0
        ( s/2, -s/2, -s/2),  # 1
        ( s/2,  s/2, -s/2),  # 2
        (-s/2,  s/2, -s/2),  # 3
        (-s/2, -s/2,  s/2),  # 4
        ( s/2, -s/2,  s/2),  # 5
        ( s/2,  s/2,  s/2),  # 6
        (-s/2,  s/2,  s/2),  # 7
    ]

    edges=[(0,1),(1,2),(2,3),(3,0),
           (4,5),(5,6),(6,7),(7,4),
           (0,4),(1,5),(2,6),(3,7)]

    points = vtk.vtkPoints()

    for i,XYZ in enumerate(local_points):
        xyz=vector['double']((XYZ[0],XYZ[1],XYZ[2]))
        p=cppyy.gbl.SPPARKS_NS.quaternion.rotate_vector(q,xyz);
        points.InsertNextPoint(x0+p[0],y0+p[1],z0+p[2])

    
    lines = vtk.vtkCellArray()
    colors = vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    namedColors=vtkNamedColors()
    for i,e in enumerate(edges):
        lines.InsertNextCell(2)
        lines.InsertCellPoint(e[0])
        lines.InsertCellPoint(e[1])
        colors.InsertNextTypedTuple(namedColors.GetColor3ub("Red"))

    # Edge ids for origin coodinate axes
    # Lower left corner
    x_id=0
    y_id=3
    z_id=8
    _b=namedColors.GetColor3ub("Blue")
    _p=namedColors.GetColor3ub("Purple")
    _g=namedColors.GetColor3ub("Green")
    colors.SetTypedTuple(x_id,_b)
    colors.SetTypedTuple(y_id,_p)
    colors.SetTypedTuple(z_id,_g)

    # Create a polydata object and set the points and lines
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetCellData().SetScalars(colors)

    # Create a mapper to map the polydata to graphics primitives
    mapper = vtk.vtkPolyDataMapper()
    # mapper.SetInputData(polydata)
    tube_filter=get_tube_filter(polydata)
    mapper.SetInputConnection(tube_filter.GetOutputPort())

    # Create an actor to represent the cube in the scene
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # Thicker lines for visibility
    actor.GetProperty().SetLineWidth(4.0)
    return actor,polydata

def create_dashed_line(p1, p2, dash_length=0.1, gap_length=0.05):
    """
    Creates a vtkPolyData representing a dashed line between two points.

    Args:
        p1 (list): Start point [x, y, z].
        p2 (list): End point [x, y, z].
        dash_length (float): Length of a visible dash segment.
        gap_length (float): Length of an invisible gap segment.

    Returns:
        vtkPolyData: PolyData containing the dashed line.
    """
    line_vector = np.array(p2) - np.array(p1)
    total_length = np.linalg.norm(line_vector)

    if total_length == 0:
        return vtk.vtkPolyData()

    unit_vector = line_vector / total_length
    current_pos = np.array(p1)

    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    segment_length = dash_length + gap_length
    num_segments = int(total_length / segment_length)
    remaining_length = total_length % segment_length

    for i in range(num_segments):
        start_dash = current_pos
        end_dash = current_pos + unit_vector * dash_length

        # Ensure we don't go past p2
        if np.linalg.norm(end_dash - np.array(p1)) > total_length:
            end_dash = np.array(p2)

        point_id1 = points.InsertNextPoint(start_dash.tolist())
        point_id2 = points.InsertNextPoint(end_dash.tolist())

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, point_id1)
        line.GetPointIds().SetId(1, point_id2)
        lines.InsertNextCell(line)

        current_pos = current_pos + unit_vector * segment_length

    # Handle any remaining part of the line that might be a dash
    if remaining_length > 0:
        start_dash = current_pos
        end_dash = current_pos + unit_vector * min(remaining_length, dash_length)

        point_id1 = points.InsertNextPoint(start_dash.tolist())
        point_id2 = points.InsertNextPoint(end_dash.tolist())

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, point_id1)
        line.GetPointIds().SetId(1, point_id2)
        lines.InsertNextCell(line)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    return polydata


def get_diagonal_actor(p0,p1):
    polydata=create_dashed_line(p0,p1)
    
    # Create a mapper to map the polydata to graphics primitives
    mapper = vtk.vtkPolyDataMapper()
    tube_filter=get_tube_filter(polydata)
    mapper.SetInputData(polydata)
    # mapper.SetInputConnection(tube_filter.GetOutputPort())

    # Create an actor to represent the cube in the scene
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # Thicker lines for visibility
    actor.GetProperty().SetLineWidth(4.0)
    actor.GetProperty().SetColor(0.0, 0.0, 0.0)
    return actor


def save_camera_state(camera):
    """Saves the essential camera parameters."""
    state = {
        "position": camera.GetPosition(),
        "focal_point": camera.GetFocalPoint(),
        "view_up": camera.GetViewUp(),
        "view_angle": camera.GetViewAngle(),
        "parallel_projection": camera.GetParallelProjection(),
        "parallel_scale": camera.GetParallelScale() if
        camera.GetParallelProjection() else None,
        "clipping_range": camera.GetClippingRange(),
        "roll": camera.GetRoll()
    }
    print("\n--- Camera State Saved ---")
    for key, value in state.items():
        print(f"{key}: {value}")
    return state

def restore_camera_state(camera, state):
    """Restores the camera to a previously saved state."""
    print("\n--- Restoring Camera State ---")
    camera.SetPosition(state["position"])
    camera.SetFocalPoint(state["focal_point"])
    camera.SetViewUp(state["view_up"])
    camera.SetViewAngle(state["view_angle"])
    camera.SetParallelProjection(state["parallel_projection"])
    if state["parallel_projection"] and state["parallel_scale"] is not None:
        camera.SetParallelScale(state["parallel_scale"])
    camera.SetClippingRange(state["clipping_range"])
    camera.SetRoll(state["roll"])
    print("Camera state restored.")

def write_camera_state(state,variant):
    fname=variant+".camera.json"
    try:
        with open(fname, 'w') as f:
            json.dump(state,f, indent=4)
            print(f"Saved current camera state to {fname}")
    except IOError as e:
        print(f"Error saving file: {e}")
    pass

def read_camera_state(variant):
    fname=variant+".camera.json"
    restored_state = {} # Initialize an empty dictionary
    try:
        with open(fname, 'r') as f:
            restored_state = json.load(f)
        print(f"\nSuccessfully restored state from {fname}")
        print(restored_state)
        print(f"Type of restored_state: {type(restored_state)}")
    except FileNotFoundError:
        print(f"Error: The file {fname} was not found.")
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from file: {e}")
    except IOError as e:
        print(f"Error reading file: {e}")
    return restored_state

def snapshot(variant,renderer):
    """
    Takes a PNG snapshot of the current rendering in a vtkRenderWindow.

    Args:
        variant (str): One of [V1, V2,...V12]; save PNG image (e.g., "V1.png").
    """
    fname=variant+".png"
    render_window=renderer.GetRenderWindow()
    
    window_to_image_filter = vtk.vtkWindowToImageFilter()
    window_to_image_filter.SetInput(render_window)
    # Capture color and alpha
    window_to_image_filter.SetInputBufferTypeToRGBA()  
    # Read from the back buffer for a clean image
    window_to_image_filter.ReadFrontBufferOff()  
    window_to_image_filter.Update()

    png_writer = vtk.vtkPNGWriter()
    png_writer.SetFileName(fname)
    png_writer.SetInputConnection(window_to_image_filter.GetOutputPort())
    png_writer.Write()

    print(f"Snapshot saved to: {fname}")

class CameraStateManager:
    def __init__(self,renderer,variant_str):
        self.variant_str=variant_str

        # Create a renderer and a render window
        self.renderer = renderer;
        self.render_window = renderer.GetRenderWindow()
        self.render_window.SetWindowName("Camera State Saver/Restorer")

        # Create an interactor
        self.render_window_interactor = vtk.vtkRenderWindowInteractor()
        self.render_window_interactor.SetRenderWindow(self.render_window)

        # # Start the interaction and rendering render_window.Render()
        style=vtk.vtkInteractorStyleTrackballCamera()
        self.render_window_interactor.SetInteractorStyle(style)

        # Camera
        self.camera=self.renderer.GetActiveCamera()
        # TODO: consider passing args to initialize camera.
        # Assume that camera is already initializeed.
        # self.renderer.ResetCamera()

        # Add key press event observer
        self.render_window_interactor.AddObserver("KeyPressEvent", 
                                                  self.key_press_callback)
        #
        print("Interact with the scene (rotate, pan, zoom).")
        print("Press 's' to save current camera state & snapshot rendering to png")
        print("Press 'r' to restore the last saved camera state.")
        print("Press 'c' to print current camera parameters.")

    def key_press_callback(self, obj, event):
        print("key_press_callback")
        key = obj.GetKeySym()
        camera = self.renderer.GetActiveCamera()

        if key == 's':
            self.save_camera_state()
            print("Camera state saved.")
        elif key == 'r':
            self.state=read_camera_state(self.variant_str)
            if self.state:
                # Restore the last saved state
                self.restore_camera_state()
                print(f"Camera state restored.")
            else:
                print("No camera state to restore. Use 's' save first.")
        elif key == 'c':
            print("\nCurrent Camera Parameters:")
            print(f"  Position: {camera.GetPosition()}")
            print(f"  Focal Point: {camera.GetFocalPoint()}")
            print(f"  View Up: {camera.GetViewUp()}")
            print(f"  Clipping Range: {camera.GetClippingRange()}")
            print(f"  Parallel Projection: {camera.GetParallelProjection()}")
            print(f"  Parallel Scale: {camera.GetParallelScale()}")
            print(f"  View Angle: {camera.GetViewAngle()}")
            print("-" * 30)
        else:
            pass

    def save_camera_state(self):
        """
        Saves the current state of a vtkCamera object into a dictionary.
        """
        camera = self.renderer.GetActiveCamera()
        self.state = {
            'Position': camera.GetPosition(),
            'FocalPoint': camera.GetFocalPoint(),
            'ViewUp': camera.GetViewUp(),
            'ClippingRange': camera.GetClippingRange(),
            'ParallelProjection': camera.GetParallelProjection(),
            'ParallelScale': camera.GetParallelScale(),
            'ViewAngle': camera.GetViewAngle()
        }
        write_camera_state(self.state,self.variant_str)
        snapshot(self.variant_str,self.renderer)
        

    def restore_camera_state(self):
        """
        Restores the state of a vtkCamera object from a dictionary.
        """
        state=self.state
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(state['Position'])
        camera.SetFocalPoint(state['FocalPoint'])
        camera.SetViewUp(state['ViewUp'])
        camera.SetClippingRange(state['ClippingRange'])
        camera.SetParallelProjection(state['ParallelProjection'])
        camera.SetParallelScale(state['ParallelScale'])
        camera.SetViewAngle(state['ViewAngle'])
        # Important to re-render after changing camera properties
        self.renderer.Render() 


    def start(self):
        self.render_window.Render()
        self.render_window_interactor.Start()

def parse_runtime_case():
    s="Render hexagonal prism orientation relative to fixed "
    s+="but oriented cube using BOR variant."
    parser = argparse.ArgumentParser(description=s)
    parser.add_argument("-r","--random", action='store_true',
                        help="Randomly orient cube")
    parser.add_argument(
        '-v', '--variant',
        choices=VALID_VARIANTS,
        required=True,
        help=f"Render specific variant. {', '.join(VALID_VARIANTS)}"
    )
    args = parser.parse_args()
    return args


# python orient_hcp.py -v 'V1'
# python orient_hcp.py -v 'V2'
# python orient_hcp.py -v 'V12'
# Randomly orient cube. (-r)
# python orient_hcp.py -v 'V2' -r
# python orient_hcp.py -v 'V12' -r
# python orient_hcp.py -v 'V9' -r
if __name__ == "__main__":
    """Burgers orientation relations for CUBIC to HCP"""
    # Parses command line
    args=parse_runtime_case()
    print(args)
    case=args.variant

    # Initial orientation of cube
    qcube=vector['double']((1,0,0,0))
    if args.random:
        # Random orienation of cube
        qcube=quaternion.generate_random_unit_quaternions(1);

    # Equivalent methods; pure python or c++ via cppyy
    variants=bor.get_variants_dictionary(qcube)
    # variants=bor.cpp_get_variants_dictionary(qcube)
    q=variants[case]

    # Point id pairs are specific to the CUBE created 
    # below. Each pair forms a diagonal line segment 
    # in the basal plane of the hexagonal prism for each 
    # of the 12 variants.  The order of these pairs are 
    # according to the ordering 'create_variant_quaternions()'

    # cases=['V'+str(i) for i in range(1,13)]
    point_pairs=[(1,7),(3,5),(0,6),(2,4),
                 (0,6),(1,7),(2,4),(3,5),
                 (3,5),(0,6),(1,7),(2,4)]
    diagonals={}
    for i,p in enumerate(point_pairs):
        diagonals['V'+str(i+1)]=p
    pair=diagonals[case]
    p0=pair[0]
    p1=pair[1]


    # z is initial orientation of 'basal' plane on HEXAGONAL prism
    z=vector['double']((0,0,1))
    # n is normal to basal plane 
    n=quaternion.rotate_vector(q,z)
    # print("%8.6f, %8.6f, %8.6f, %8.6f"%(q[0], q[1], q[2], q[3]))

    # dimension
    h=1.0
    # position
    c=(0,0,0)

    # HEXAGONAL prism
    side,height=h*sqrt(3),2.0*h
    args={'orient':q,'center':c,'side':side,'height':height}
    hex_actor=get_hexagonal_prism_actor(**args)

    # Move HEX to align points with fixed but oriented cube.
    hex_actor.AddPosition(tuple(-x*h for x in n))

    # Oriented CUBE
    args={'orient':qcube,'center':c,'edge':2*h}
    cube_actor,cube_polydata=get_cube_actor(**args)
    cube_mapper=cube_actor.GetMapper()
    # cube_polydata=cube_mapper.GetInput()
    # cube_points=cube_polydata.GetPoints()
    # a=cube_points.GetPoint(p0)
    # b=cube_points.GetPoint(p1)
    # diag_actor=get_diagonal_actor(a,b)

    # Create basal plane, using normal
    plane_actor=get_plane_actor(c,(n[0],n[1],n[2]),size=h)

    # Create renderer and render window
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetWindowName("HCP orientation wrt CUBIC")
    desired_width=1200
    desired_height=800
    render_window.SetSize(desired_width,desired_height)

    # Add the actors to the renderer
    # renderer.AddActor(plane_actor)
    renderer.AddActor(cube_actor)
    renderer.AddActor(hex_actor)
    # renderer.AddActor(diag_actor)

    # Set background color
    renderer.SetBackground(1.0, 1.0, 1.0)

    # Initial camera orientation etc.
    camera=renderer.GetActiveCamera()
    camera.SetPosition(tuple(15*x for x in n))
    camera.SetFocalPoint(tuple(0*x for x in n))
    camera.SetViewUp(0.0, 0.0, 1.0)

    state_manager=CameraStateManager(renderer,case)
    state_manager.start()
