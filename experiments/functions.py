#############################################################################################################################
######################################################## Global stuff #######################################################
#############################################################################################################################

# Imports
import pandas
import cairo
import igraph
import itertools
import scipy.spatial as spatial
import matplotlib.pyplot as pyplot
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import numpy
import numpy.linalg as linalg
import ot
import os
import IPython.display as display
import scipy.cluster.hierarchy as hierarchy
import sklearn.cluster as skcluster
import sklearn.manifold as manifold
import scipy.spatial.distance as ssd
import lifelines
import lifelines.statistics as statistics
import hashlib
import sys
import scipy.stats as stats
import pickle
import heapq
import copy
import multiprocessing
import inspect
import cProfile
import statistics as sts 



# Constants
CATEGORIES = {1 : "Lymphocyte", 2 : "Stroma", 3 : "Cancer"}
RATIO_PX_TO_UM = {40 : 1262.5/5000, 20 : 1481.7/3000} 
MEMO_DIR = "memo"
DEFAULT_BINS = 100
NB_PROCESSES = 8

#############################################################################################################################
###################################################### Utlity functions #####################################################
#############################################################################################################################

"""
    Starts profiling execution times in the code.
"""

def start_profiling () :
    
    # Create global object
    global profiler
    profiler = cProfile.Profile()
    profiler.enable()

#############################################################################################################################

"""
    Stops profiling execution times in the code, and shows results.
"""

def stop_profiling (file_name=None) :
    
    # Get stats
    global profiler
    profiler.create_stats()
    times = []
    functions = []
    for function_called in profiler.getstats() :
        if "functions.py" in str(function_called.code) :
            function_label = str(function_called.code).split("<code ")[1].split(" at ")[0] + " (" + str(function_called.callcount) + ")"
            functions.append(function_label)
            times.append(function_called.totaltime)
            
    # Plot
    figure = pyplot.figure(figsize=(20, 10))
    pyplot.bar(range(len(functions)), times)
    pyplot.xticks(range(len(functions)), functions, rotation=90)
    pyplot.show()

    # Save
    if file_name is not None :
        create_directory_for(file_name)
        figure.savefig(file_name)

#############################################################################################################################

"""
    Loads the memoized result of a function.
"""

def memoize_load (*args) :
    
    # Save in a directory corresponding with the function call
    caller_function = sys._getframe().f_back.f_code.co_name
    memo_dir = MEMO_DIR + "/" + caller_function
    
    # Create hash string of the id string "str(arg1)_str(arg2)_..._str(argX)"
    # For lambda functions we use the code
    id_string = "_".join([inspect.getsource(arg) if callable(arg) else str(arg) for arg in args])
    hash_string = hashlib.md5(id_string.encode()).hexdigest()
    memo_file_name = memo_dir + "/" + hash_string + ".pkl"

    # Return existing file contents or False
    create_directory_for(memo_file_name)
    if os.path.exists(memo_file_name) :
        with open(memo_file_name, "rb") as memo_file:
            contents = pickle.load(memo_file)
        return True, contents
    return False, memo_file_name

#############################################################################################################################

"""
    Memoizes the result of a function.
"""

def memoize_save (memo_file_name, variable_to_save) :
    
    # Save args if asked
    with open(memo_file_name, "wb") as memo_file :
        pickle.dump(variable_to_save, memo_file)
    
#############################################################################################################################

"""
    Creates the directory for the given file name if it does not exist.
"""

def create_directory_for (file_name) :
    
    # Creates the corresponding directory
    dir_name = os.path.dirname(file_name)
    os.makedirs(dir_name, exist_ok=True)
    
#############################################################################################################################

"""
    Loads a ROI from a given file.
    We need to precise the magnification to work in µm and not in px.
"""

def load_roi (file_name, magnification=40) :
    
    # Load first sheet from the file
    file_object = pandas.ExcelFile(file_name)
    data_frame = file_object.parse(file_object.sheet_names[0])

    # Update cell coordinates to work in µm
    data_frame["x_um"] = data_frame["Location.X"] * RATIO_PX_TO_UM[magnification]
    data_frame["y_um"] = data_frame["Location.Y"] * RATIO_PX_TO_UM[magnification]
    data_frame = data_frame.rename(columns={"Location.X" : "x_px", "Location.Y" : "y_px"})
    return data_frame

#############################################################################################################################

"""
    Loads the clinical information for the given patients.
"""

def load_clinical (patient_ids, file_name) :
    
    # Load first sheet from the file and keep the good rows
    data_frame = pandas.read_csv(file_name)
    data_frame = data_frame[data_frame["ID"].isin(patient_ids)]
    data_frame = data_frame.replace({"Vital_Status": {"Alive": 0, "Dead": 1}})
    return data_frame


#############################################################################################################################

"""
    Condition to select valid ROIs : 
    (0.3 * Number of cells < Number of Cancer cells < 0.7 * Number of Cells) & ( Number of Stroma cells > 0.05 * Number of cells).
    The idea is to select only ROI with multiples types of cells, so we can study their interactions.
"""

def valid_roi(cells):
    number_of_cancer_cells=cells[cells['Category']==3].count()[0]
    number_of_stroma_cells=cells[cells['Category']==2].count()[0]
    number_of_lymphocyte_cells=cells[cells['Category']==1].count()[0]
    number_of_cells=cells.count()[0]
    if (30*number_of_cells/100)<number_of_cancer_cells<(70*number_of_cells/100) and number_of_stroma_cells>(5*number_of_cells/100):
        return True
    else :
        return False

#############################################################################################################################

"""
    Builds a Delaunay triangulation for a data frame of cells.
"""

def delaunay_triangulation (data_frame, memoize=True) :
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(data_frame)
        if loaded :
            return contents_or_target
    
    # Triangulate
    points = data_frame[["x_um", "y_um"]].values.tolist()
    triangulation = spatial.Delaunay(points)
    
    # Build graph
    matrix = numpy.zeros((data_frame.shape[0], data_frame.shape[0]))
    for i in range(data_frame.shape[0]):
        for j in triangulation.vertex_neighbor_vertices[1][triangulation.vertex_neighbor_vertices[0][i] : triangulation.vertex_neighbor_vertices[0][i + 1]] :
            matrix[i, j] = 1
    graph = igraph.Graph.Adjacency((matrix > 0).tolist(), igraph.ADJ_UNDIRECTED)
    
    # We save the data frame and triangulation in the graph object
    graph.data_frame = data_frame
    graph.triangulation = triangulation
    
    # Append Euclidean weights to edges
    points_1 = [edge.tuple[0] for edge in graph.es]
    points_2 = [edge.tuple[1] for edge in graph.es]
    distances = euclidean_distance(graph.data_frame, points_1, points_2)
    for i in range(len(graph.es)) :
        graph.es[i]["weight"] = distances[i]

    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, graph)
    return graph
    
#############################################################################################################################

"""
    Finds the triangles within the Delaunay triangulation that interface particular vertices.
"""

def find_interface_triangles (graph, is_vertex_of_type_1, is_vertex_of_type_2, is_vertex_of_type_3=None) :
    
    # We iterate over all triangles to find the ones that link vertices of the 3 specified types
    nb_types = 2 if is_vertex_of_type_3 is None else 3
    vertex_types = {vertex : [is_vertex_of_type_1(vertex), is_vertex_of_type_2(vertex), True if nb_types == 2 else is_vertex_of_type_3(vertex)] for vertex in range(graph.vcount())}
    def check_triangle (triangle_index) :
        triangle = graph.triangulation.simplices[triangle_index]
        if all(vertex < graph.vcount() for vertex in triangle) :
            for vertices in itertools.permutations(triangle, nb_types) :
                if vertex_types[vertices[0]][0] and vertex_types[vertices[1]][1] :
                    if nb_types == 2 or (nb_types == 3 and vertex_types[vertices[2]][2]) :
                        return triangle_index
        return None
    with multiprocessing.pool.ThreadPool(processes=NB_PROCESSES) as pool :
        results = pool.map(check_triangle, range(len(graph.triangulation.simplices)))
        interface_triangles = [result for result in results if result is not None]
        triangle_index_to_list_index = {interface_triangles[i] : i for i in range(len(interface_triangles))}
    
    # We create a graph of triangles (there is an edge if they share 2 vertices)
    triangles_graph = igraph.Graph()
    triangles_graph.add_vertices(len(interface_triangles))
    for list_i in range(len(interface_triangles)) :
        triangle_i = interface_triangles[list_i]
        for triangle_j in graph.triangulation.neighbors[triangle_i] :
            if triangle_i > triangle_j :
                if triangle_j in triangle_index_to_list_index :
                    list_j = triangle_index_to_list_index[triangle_j]
                    triangles_graph.add_edges([(list_i,list_j)])
    
    # We group interface triangles by connected component in that graph
    connected_components = triangles_graph.clusters()
    clusters_of_triangles = []
    for component in connected_components :
        component_triangles = []
        for triangle_id in component :
            component_triangles.append(list(graph.triangulation.simplices[interface_triangles[triangle_id]]))
        clusters_of_triangles.append(component_triangles)
    return clusters_of_triangles

#############################################################################################################################

"""
    Finds the groups of a particular type of cells.
"""

def find_cell_type_groups (graph, is_vertex_of_type) :
    
    # We keep edges that link pairs of vertices of that type
    edge_type_graph = igraph.Graph()
    edge_type_graph.add_vertices(graph.vcount())
    for edge in graph.es :
        if all(is_vertex_of_type(vertex) for vertex in edge.tuple) :
            edge_type_graph.add_edges([edge.tuple])
    
    # We find the remaining connected components (and keep those of the vertex type of interest)
    connected_components = edge_type_graph.clusters()
    clusters_of_cell_type = []
    for component in connected_components :
        if len(component) > 1 or (len(component) == 1 and is_vertex_of_type(component[0])) :
            clusters_of_cell_type.append(component)
    return clusters_of_cell_type

#############################################################################################################################

"""
    Computes the area of a given shape described by points.
"""

def shape_area (points, coordinates_dataframe) :
    
    # Shoelace formula
    x = coordinates_dataframe.iloc[points]["x_um"]
    y = coordinates_dataframe.iloc[points]["y_um"]
    return 0.5 * numpy.abs(numpy.dot(x, numpy.roll(y, 1)) - numpy.dot(y, numpy.roll(x, 1)))

#############################################################################################################################

"""
    Plots a graph.
"""

def plot_graph (graph, signal=None, legend=None, shapes=None, file_name=None) :
    
    # Create output directory if needed
    if file_name is None :
        file_name = "tmp.png"
    else :
        create_directory_for(file_name)
    
    # Global visual attributes
    visual_style = {}
    visual_style["vertex_size"] = 10
    visual_style["bbox"] = (3000, 1500)
    visual_style["margin"] = 10
    visual_style["background"] = [0.0, 0.0, 0.0, 0.0]
    visual_style["edge_curved"] = False
    
    # We need coordinates
    graph_to_plot = copy.copy(graph)
    for vertex in range(graph_to_plot.vcount()) :
        graph_to_plot.vs[vertex]["x"] = graph_to_plot.data_frame.iloc[vertex]["x_px"]
        graph_to_plot.vs[vertex]["y"] = graph_to_plot.data_frame.iloc[vertex]["y_px"]
    
    # Add colors if asked
    if signal is None :
        for vertex in range(graph_to_plot.vcount()) :
            graph_to_plot.vs[vertex]["color"] = [1.0, 1.0, 1.0]
    elif legend is None :
        if any(value < 0.0 or value > 1.0 for value in signal) :
            print("Warning: signal contains entries not in [0, 1], normalizing (-min then /max)")
            signal = [(value - min(signal)) / max(signal) for value in signal]
        for vertex in range(graph_to_plot.vcount()) :
            graph_to_plot.vs[vertex]["color"] = [signal[vertex], 0.5, 1.0 - signal[vertex]]
    else :
        palette = igraph.ClusterColoringPalette(len(legend))
        sorted_legend = sorted(legend.values())
        colormap = {value : palette[sorted_legend.index(legend[value])] for value in legend}
        for vertex in range(graph_to_plot.vcount()) :
            graph_to_plot.vs[vertex]["color"] = colormap[signal[vertex]]
    
    # Export
    igraph.plot(graph_to_plot, file_name, **visual_style)
    
    # Add legend if asked
    if legend is not None :
        font_size = 20
        surface = cairo.ImageSurface.create_from_png(file_name)
        context = cairo.Context(surface)
        context.set_font_size(font_size)
        max_text_width = max([context.text_extents(value).width for value in sorted_legend])
        context.set_line_width(2)
        context.set_source_rgba(1.0, 1.0, 1.0, 0.8)
        context.rectangle(font_size, font_size, font_size + max_text_width, font_size * (0.5 + len(sorted_legend)))
        context.fill()
        context.set_source_rgba(0.0, 0.0, 0.0, 0.8)
        context.rectangle(font_size, font_size, font_size + max_text_width, font_size * (0.5 + len(sorted_legend)))
        context.stroke()
        y = 2.0 * font_size
        for i in range(len(sorted_legend)) :
            context.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            for x_shift in [-1, 0, 1] : 
                for y_shift in [-1, 0, 1] : 
                    context.move_to(1.5 * font_size + x_shift, y + y_shift)
                    context.show_text(sorted_legend[i])
            context.set_source_rgba(*colormap[list(legend.keys())[list(legend.values()).index(sorted_legend[i])]])
            context.move_to(1.5 * font_size, y)
            context.show_text(sorted_legend[i])
            y += font_size
        context.stroke()
        surface.write_to_png(file_name)
    
    # Add surfaces if asked
    if shapes is not None :
        min_x = min([graph_to_plot.vs[i]["x"] for i in range(len(graph_to_plot.vs))])
        max_x = max([graph_to_plot.vs[i]["x"] for i in range(len(graph_to_plot.vs))])
        min_y = min([graph_to_plot.vs[i]["y"] for i in range(len(graph_to_plot.vs))])
        max_y = max([graph_to_plot.vs[i]["y"] for i in range(len(graph_to_plot.vs))])
        nb_clusters = len(shapes)
        shapes_palette = igraph.ClusterColoringPalette(nb_clusters)
        surface = cairo.ImageSurface.create_from_png(file_name)
        surface_back = cairo.ImageSurface(surface.get_format(), surface.get_width(), surface.get_height())
        context = cairo.Context(surface_back)
        context.set_line_width(6)
        for cluster_id in range(nb_clusters) :
            cluster = shapes[cluster_id]
            context.set_source_rgba(*shapes_palette[cluster_id])
            for shape in cluster :
                for point in shape :
                    point_x = (graph_to_plot.vs[point]["x"] - min_x) * (visual_style["bbox"][0] - 2 * visual_style["margin"]) / (max_x - min_x) + visual_style["margin"]
                    point_y = (graph_to_plot.vs[point]["y"] - min_y) * (visual_style["bbox"][1] - 2 * visual_style["margin"]) / (max_y - min_y) + visual_style["margin"]
                    if point == shape[0] :
                        context.move_to(point_x, point_y)
                    else :
                        context.line_to(point_x, point_y)
                if len(shape) > 2 :
                    context.close_path()
                    context.fill()
                context.stroke()
        context.set_source_surface(surface)
        context.paint()
        surface_back.write_to_png(file_name)
    
    # Plot
    pyplot.figure(figsize=(20, 10))
    pyplot.axis("off")
    pyplot.imshow(pyplot.imread(file_name))
    pyplot.tight_layout()
    pyplot.show()
    
    # Remove tmp file if there is one
    if os.path.exists("tmp.png") :
        os.remove("tmp.png")
    
#############################################################################################################################

"""
    Dijkstra's algorithm that stops when a condition is verified.
    This returns the closest vertex according to the Euclidean distance.
    This work only because Delaunay triangulation is based on Euclidean coordinates and triangular inequality is verified.
"""

def dijkstra_find_closest (graph, start_vertex, condition) :
    
    # Initialize with starting vertex
    queuing_structure = []
    heapq.heappush(queuing_structure, (0, start_vertex))
    explored_vertices = {}
    
    # Traversal
    graph_distance_threshold = None
    kept_vertices = []
    while len(queuing_structure) > 0 :
        
        # Pop element
        (graph_distance_to_current_vertex, current_vertex) = heapq.heappop(queuing_structure)
        
        # If we have went far enough, we return the closest vertex that verified the condition (using Euclidean distance)
        if graph_distance_threshold is not None and graph_distance_to_current_vertex > graph_distance_threshold :
            return heapq.heappop(kept_vertices)
        
        # If this validates the condition, we retain it
        # If it does so for the first time, we set the stop threshold by multiplying the Euclidean distance by 4 pi / 3 sqrt(3)
        # https://en.wikipedia.org/wiki/Delaunay_triangulation
        # Keil, J. Mark; Gutwin, Carl A. (1992), "Classes of graphs which approximate the complete Euclidean graph"
        if condition(graph, start_vertex, current_vertex) :
            euclidean_vertex_distance = euclidean_distance(graph.data_frame, start_vertex, current_vertex)
            heapq.heappush(kept_vertices, (euclidean_vertex_distance, current_vertex))
            if graph_distance_threshold is None :
                graph_distance_threshold = euclidean_vertex_distance * 4 * numpy.pi / (3 * numpy.sqrt(3))
        
        # Consider neighbors
        if current_vertex not in explored_vertices :
            explored_vertices[current_vertex] = graph_distance_to_current_vertex
            for neighbor in graph.neighbors(current_vertex) :
                if neighbor not in explored_vertices :
                    edge_weight = graph.es.select(_source=current_vertex, _target=neighbor)["weight"][0]
                    graph_distance_to_neighbor = graph_distance_to_current_vertex + edge_weight
                    heapq.heappush(queuing_structure, (graph_distance_to_neighbor, neighbor))
    
    # Default return if not found
    return float("inf"), None

#############################################################################################################################

"""
    Dijkstra's algorithm that stops at some distance.
    This returns all vertices below that distance considering Euclidean distance.
    Uses the same property as above with the Delaunay triangulation.
"""

def dijkstra_find_all_vertices_closer_than (graph, start_vertex, condition, euclidean_threshold) :
    
    # Initialize with starting vertex
    queuing_structure = []
    heapq.heappush(queuing_structure, (0, start_vertex))
    explored_vertices = {}
    
    # Traversal
    graph_distance_threshold = euclidean_threshold * 4 * numpy.pi / (3 * numpy.sqrt(3))
    kept_vertices = []
    while len(queuing_structure) > 0 :
        
        # Pop element
        (graph_distance_to_current_vertex, current_vertex) = heapq.heappop(queuing_structure)
        
        # If we have went far enough, we return the closest vertex that verified the condition (using Euclidean distance)
        if graph_distance_to_current_vertex > graph_distance_threshold :
            return kept_vertices
        
        # If this validates the condition, we retain it
        if condition(graph, start_vertex, current_vertex) :
            euclidean_vertex_distance = euclidean_distance(graph.data_frame, start_vertex, current_vertex)
            heapq.heappush(kept_vertices, (euclidean_vertex_distance, current_vertex))
        
        # Consider neighbors
        if current_vertex not in explored_vertices :
            explored_vertices[current_vertex] = graph_distance_to_current_vertex
            for neighbor in graph.neighbors(current_vertex) :
                if neighbor not in explored_vertices :
                    edge_weight = graph.es.select(_source=current_vertex, _target=neighbor)["weight"][0]
                    graph_distance_to_neighbor = graph_distance_to_current_vertex + edge_weight
                    heapq.heappush(queuing_structure, (graph_distance_to_neighbor, neighbor))
    
    # Default return if not found
    return []

#############################################################################################################################

"""
    Euclidean distances between vertices_1[i] and vertices_2[i].
"""

def euclidean_distance (data_frame, vertices_1, vertices_2) :
    
    # L2 norm of difference
    coordinates_1 = numpy.concatenate((data_frame.iloc[vertices_1]["x_um"], data_frame.iloc[vertices_1]["y_um"])).reshape(2, -1)
    coordinates_2 = numpy.concatenate((data_frame.iloc[vertices_2]["x_um"], data_frame.iloc[vertices_2]["y_um"])).reshape(2, -1)
    difference = coordinates_1 - coordinates_2
    return linalg.norm(difference, axis=0)

################################################################################################################

"""
    Euclidean distances between all pairs of vertices (v1 in vertices_1, v2 in vertices_2).
"""

def all_pairwise_euclidean_distances (data_frame, vertices_1, vertices_2) :
    
    # L2 norm of difference
    coordinates_1 = numpy.concatenate((data_frame.iloc[vertices_1]["x_um"], data_frame.iloc[vertices_1]["y_um"])).reshape(2, -1)
    coordinates_2 = numpy.concatenate((data_frame.iloc[vertices_2]["x_um"], data_frame.iloc[vertices_2]["y_um"])).reshape(2, -1)
    difference = coordinates_1[:, numpy.newaxis, :] - coordinates_2[:, :, numpy.newaxis]
    return linalg.norm(difference, axis=0).transpose()

#############################################################################################################################

"""
    Returns shortest distances between some given types of vertices, from type 1 to type 2.
"""

def shortest_distances_between_types (graph, is_vertex_of_type_1, is_vertex_of_type_2) :

    # Identify vertices of types
    type_1_vertices = [vertex for vertex in range(graph.vcount()) if is_vertex_of_type_1(vertex)]
    if len(type_1_vertices) == 0 :
        return []
    type_2_vertices = [vertex for vertex in range(graph.vcount()) if is_vertex_of_type_2(vertex)]
    if len(type_2_vertices) == 0 :
        return [float("inf")] * len(type_1_vertices)
    
    # Compute all pairwise distances and return min per row
    all_distances = all_pairwise_euclidean_distances(graph.data_frame, type_1_vertices, type_2_vertices)
    return numpy.min(all_distances, axis=1).tolist()

#############################################################################################################################

"""
    Transforms a list of values into a normalized histogram.
"""

def to_distribution (values, min_x=None, max_x=None, nb_bins=DEFAULT_BINS) :

    # Compute histogram
    if min_x is None :
        min_x = min(values)
    if max_x is None :
        max_x = max(values)
    all_values_hist, all_values_bins = numpy.histogram(values, bins=nb_bins, range=(min_x, max_x))
    
    # Normalize
    ratio_hists = [all_values_hist[i] / len(values) for i in range(len(all_values_hist))]
    return all_values_bins, ratio_hists

#############################################################################################################################

"""
    Plots a distribution.
"""

def plot_distribution (values, min_x=None, max_x=None, nb_bins=DEFAULT_BINS, xlabel="", file_name=None) :
        
    # Compute distribution
    if min_x is None :
        min_x = min(values)
    if max_x is None :
        max_x = max(values)
    bins, probas = to_distribution(values, min_x, max_x, nb_bins=nb_bins)
    
    # Plot histogram
    figure = pyplot.figure(figsize=(20, 10))
    pyplot.bar(bins[:-1], probas, width=(bins[1] - bins[0]), align="edge")
    pyplot.xlabel(xlabel)
    pyplot.ylabel("%")
    pyplot.ylim([0, 1])
    pyplot.xlim([min_x, max_x])
    pyplot.tight_layout()
    pyplot.show()
    
    # Save
    if file_name is not None :
        create_directory_for(file_name)
        figure.savefig(file_name)

#############################################################################################################################

"""
    Plots a matrix.
"""

def plot_matrix (matrix, colorbar=False, round_values=None, file_name=None) :

    # Plot matrix
    figure, axis = pyplot.subplots(figsize=(20, 20))
    axis.matshow(matrix)
    
    # Add colorbar
    if colorbar :
        pyplot.colorbar()
    
    # Add values
    if round_values is not None :
        for i in range(matrix.shape[0]) :
            for j in range(matrix.shape[1]) :
                value = round(matrix[i, j], round_values)
                color = "black" if matrix[i, j] > 0.5 else "white"
                axis.text(i, j, str(value), va="center", ha="center", color=color)
    pyplot.show()
    
    # Save
    if file_name is not None :
        create_directory_for(file_name)
        mpimg.imsave(file_name, matrix)
        
#############################################################################################################################

"""
    Computes the earth mover distance between two distance profiles.
"""

def wasserstein (values_1, values_2) :
    
    # Degenerate cases
    if len(values_1) == 0 or len(values_2) == 0 :
        return float("inf")
    
    # Wasserstein distance
    return stats.wasserstein_distance(values_1, values_2)
    
#############################################################################################################################

"""
    Computes the Kolmogorov-Smirnov test between two distance profiles.
"""

def kolmogorov_smirnov (values_1, values_2) :
    
    # Degenerate cases
    if len(values_1) == 0 or len(values_2) == 0 :
        return float("inf")
    
    # Perform test
    ks_stats, p_value = stats.ks_2samp(values_1, values_2)
    return ks_stats,p_value

#############################################################################################################################

"""
    Computes the Kullback-Leibler divergence between two distance profiles.
"""

def kullback_leibler (values_1, values_2, min_x=None, max_x=None, nb_bins=DEFAULT_BINS) :
    
    # Degenerate cases
    if len(values_1) == 0 or len(values_2) == 0 :
        return float("inf")
    
    # We want to set the same min_x and max_x when comparing more than 2 histograms (to have the same bins)
    if min_x is None :
        min_x = min(values_1 + values_2)
    if max_x is None :
        max_x = max(values_1 + values_2)
    
    # Compute value
    bins, probas_1 = to_distribution(values_1, min_x, max_x, nb_bins=nb_bins)
    bins, probas_2 = to_distribution(values_2, min_x, max_x, nb_bins=nb_bins)
    probas_1 = numpy.array(probas_1)
    probas_2 = numpy.array(probas_2)
    probas_1[numpy.argwhere(probas_2 == 0)] = 0
    kl_divergence = stats.entropy(probas_1, probas_2)
    return kl_divergence

#############################################################################################################################

"""
    Flattens an arbitrarily nested containers.
"""

def flatten (container) :
    
    # If not a container, we return it
    if not isinstance(container, (list, tuple)) :
        return [container]
    
    # Otherwise we flatten
    flattened_container = []
    for sub_container in container :
        flattened_container += flatten(sub_container)
    return flattened_container

#############################################################################################################################

"""
    Performs a hierarchical clustering using a distances matrix.
"""

def hierarchical_clustering (distances, labels, nb_clusters=2, file_name=None) :
    
    # If there are infinite values, we add a cluster for them (will be cluster 0) and then remove them
    clusters = [-1] * distances.shape[0]
    outlier_indices = numpy.where(numpy.sum(numpy.isinf(distances), axis=0) >= len(labels) - 1)[0]
    ok_indices = [i for i in range(distances.shape[0]) if i not in outlier_indices]
    corrected_distances = numpy.delete(distances, outlier_indices, axis=0)
    corrected_distances = numpy.delete(corrected_distances, outlier_indices, axis=1)
    corrected_labels = numpy.delete(labels, outlier_indices)
    
    # Perform clustering
    model = skcluster.AgglomerativeClustering(affinity="precomputed", linkage="single", n_clusters=nb_clusters).fit(corrected_distances)
    for i in range(corrected_distances.shape[0]) :
        clusters[ok_indices[i]] = model.labels_[i]
    
    # Find max distance within each cluster below threshold,
    # or min non-null distance within each cluster above it
    max_distances = []
    for cluster_id in range(nb_clusters) :
        if numpy.count_nonzero(numpy.array(clusters) == cluster_id) == 1 :
            submatrix_index = int(numpy.where(numpy.array(clusters) == cluster_id)[0])
            submatrix_distances = distances[:, submatrix_index]
            submatrix_distances[submatrix_index] = -1
            submatrix_min = numpy.min(submatrix_distances[submatrix_distances >= 0])
            max_distances.append(submatrix_min)
        else :
            submatrix_indices = [i for i in range(len(clusters)) if clusters[i] == cluster_id]
            submatrix_distances = distances[numpy.ix_(submatrix_indices, submatrix_indices)]
            submatrix_max = numpy.max(submatrix_distances)
            max_distances.append(submatrix_max)
    
    # Sort clusters in decreasing max order
    max_distances_order = numpy.argsort(max_distances)
    max_distances_order = numpy.flip(max_distances_order)
    for i in range(len(clusters)) :
        if clusters[i] != -1 :
            clusters[i] = int(numpy.where(max_distances_order == clusters[i])[0])
        if corrected_distances.shape != distances.shape :
            clusters[i] += 1
    
    # Plot dendrogram if asked
    if file_name is not None :
        create_directory_for(file_name)
        linkage_matrix = hierarchy.linkage(ssd.squareform(corrected_distances), "single")
        figure = pyplot.figure(figsize=(20, 10))
        hierarchy.dendrogram(linkage_matrix, labels=corrected_labels, leaf_rotation=90, color_threshold=linkage_matrix[:, 2][-nb_clusters+1])
        pyplot.tight_layout()
        pyplot.show()
        figure.savefig(file_name)

    # Return labels in the clusters
    clusters_of_labels = [[] for i in range(len(set(clusters)))]
    for i in range(len(clusters)) :
        clusters_of_labels[clusters[i]].append(labels[i])
    return clusters_of_labels

#############################################################################################################################

"""
    Performs a spectral clustering using a distances matrix.
"""

def spectral_clustering (distances, labels, nb_clusters=2, file_name=None) :
    
    # If there are infinite values, we add a cluster for them (id nb_cluster) and then remove them
    clusters = [nb_clusters] * distances.shape[0]
    outlier_indices = numpy.where(numpy.sum(numpy.isinf(distances), axis=0) >= len(labels) - 1)[0]
    ok_indices = [i for i in range(distances.shape[0]) if i not in outlier_indices]
    corrected_distances = numpy.delete(distances, outlier_indices, axis=0)
    corrected_distances = numpy.delete(corrected_distances, outlier_indices, axis=1)
    corrected_labels = numpy.delete(labels, outlier_indices)
    
    # We work on the RBF of the distances matrix
    affinity_matrix = numpy.exp(-corrected_distances)
    
    # Perform clustering
    model = skcluster.SpectralClustering(affinity="precomputed", n_clusters=nb_clusters).fit(affinity_matrix)
    for i in range(affinity_matrix.shape[0]) :
        clusters[ok_indices[i]] = model.labels_[i]

    # Plot projection if asked
    if file_name is not None :
        create_directory_for(file_name)
        coords = manifold.SpectralEmbedding(affinity="precomputed", n_components=2).fit_transform(affinity_matrix)
        pyplot.figure(figsize=(20, 10))
        pyplot.scatter(coords[:, 0], coords[:, 1], c=numpy.array(clusters)[ok_indices])
        pyplot.tight_layout()
        pyplot.show()
        pyplot.savefig(file_name)

    # Return labels in the clusters
    clusters_of_labels = [[] for i in range(len(set(clusters)))]
    for i in range(len(clusters)) :
        clusters_of_labels[clusters[i]].append(labels[i])
    return clusters_of_labels

#############################################################################################################################

"""
    Plots the Kaplan-Meier survival curves for all given populations.
"""

def plot_survival_curves (populations, clinical_data, file_name=None) :

    # We plot all curves in a single figure
    figure = pyplot.figure(figsize=(6, 4))
    ax = None
    for i in range(len(populations)) :
        population_data = clinical_data[clinical_data["ID"].isin(populations[i])]
        fitter = lifelines.KaplanMeierFitter()
        fitter.fit(population_data["Days_to_last_follow_up"], population_data["Vital_Status"], label="Population " + str(2-i) + " (" + str(population_data.shape[0]) + ")")
        if figure :
            fitter.plot(ax=ax, show_censors=True)
        else :
            ax = fitter.plot(show_censors=True)
    
    # Additional parameters
    all_populations_data = clinical_data[clinical_data["ID"].isin(flatten(populations))]
    pyplot.xlim(min(all_populations_data["Days_to_last_follow_up"]), max(all_populations_data["Days_to_last_follow_up"]))
    pyplot.ylim(0.0, 1.0)
    pyplot.xlabel('TTLFU (days)')
    pyplot.ylabel('% alive')
    #pyplot.title("Kaplan-Meier")
    pyplot.tight_layout()
    pyplot.show()
    
    # Save
    if file_name is not None :
        create_directory_for(file_name)
        figure.savefig(file_name)
        
#############################################################################################################################

"""
    Performs log-rank test between populations (pairwise).
"""

def log_rank_test (populations, clinical_data) :

    # We perform the log-rank test
    results = numpy.zeros((len(populations), len(populations)))
    p_values = numpy.zeros((len(populations), len(populations)))
    for i in range(len(populations)) :
        data_frame_i = clinical_data[clinical_data["ID"].isin(populations[i])]
        for j in range(i + 1, len(populations)) :
            data_frame_j = clinical_data[clinical_data["ID"].isin(populations[j])]
            test_result = statistics.logrank_test(data_frame_i["Days_to_last_follow_up"], data_frame_j["Days_to_last_follow_up"], data_frame_i["Vital_Status"], data_frame_j["Vital_Status"])
            results[i, j] = test_result.test_statistic
            results[j, i] = results[i, j]
            p_values[i, j] = test_result.p_value
            p_values[j, i] = p_values[i, j]
    return results, p_values

#############################################################################################################################

"""
    Used in density feature.
"""

def Convert(lst): 
    return [ -i for i in lst ] 


#############################################################################################################################

"""
    Used in density feature. Merge dictionnaries
"""

def merge_dictionaries(dict1, dict2):
    merged_dictionary = {}
    
    for key in dict1:
        if key in dict2:
            new_value = dict1[key]+dict2[key]
        else:
            new_value = dict1[key]

        merged_dictionary[key] = new_value

    for key in dict2:
        if key not in merged_dictionary:
            merged_dictionary[key] = dict2[key]

    return merged_dictionary

#############################################################################################################################

"""
    Used in density feature. Return a dict of mean values per keys.
"""

def mean_dict(dict_):
    meaned_dictionary = {}
    for key in dict_:
        meaned_dictionary[key]=sts.mean(dict_[key])
    return meaned_dictionary

#############################################################################################################################

"""
    Compute Wasserstein distance between histograms.
"""

def wasserstein_between_histograms (bins_1, values_1, bins_2, values_2) :
    # If one set does not contain any value, we set an infinite cost
    if sum(values_1) == 0 or sum(values_2) == 0 :
        return float("inf")
    
    # Wasserstein distance
    return stats.wasserstein_distance(bins_1, bins_2, values_1, values_2)



#############################################################################################################################

"""
    Plots a curve.
"""

def plot_curve (xs, ys, xlabel="", ylabel="", file_name=None) :
        
    # Plot
    figure = pyplot.figure(figsize=(20, 10))
    pyplot.plot(xs, ys)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.tight_layout()
    pyplot.show()
    
    # Save
    if file_name is not None :
        create_directory_for(file_name)
        figure.savefig(file_name)
#############################################################################################################################

"""
    Plots central cluster patient histogram.
"""

def plot_central_histogram(feature_name,patient_id,num_cluster):
    if feature_name !='lymphocyte_density':
        pyplot.figure(figsize=(7.5,4.8))
        pyplot.hist(all_features[patient_id][feature_name],len(all_features[patient_id][feature_name]))
        if num_cluster == '0':
            num = '2'
        else :
            num = '1'
        pyplot.xlabel(feature_name+ ' S'+num)
        pyplot.ylabel('Number of occurrences')
        pyplot.show()
        pyplot.savefig(feature_name+'_cluster'+num+'.png')
        
#############################################################################################################################

"""
    Compute diverse statistics about the obtained clusters.
"""       

def compute_cluster_statistics(feature_name,type_):

    list_moy0_f=[]
    list_moy1_f=[]
    if type_== 'small':
        var = 0
    if type_ == 'big':
        var = -1
    if feature_name != 'lymphocyte_density':

        for patient_id in dict_clusters[feature_name][0]:
            all_features[patient_id][feature_name].sort()
            if all_features[patient_id][feature_name] != []:
                list_moy0_f.append(all_features[patient_id][feature_name][var])

        list_moy0_f = [k for k in list_moy0_f if k != numpy.inf]


        for patient_id in dict_clusters[feature_name][1]:
            all_features[patient_id][feature_name].sort()
            if all_features[patient_id][feature_name] != []:
                list_moy1_f.append(all_features[patient_id][feature_name][var])

        list_moy1_f = [k for k in list_moy1_f if k != numpy.inf]
        
        print('Mean '+type_+'est distance S2 : ',sts.mean(list_moy0_f))
        print('Mean '+type_+'est distance S1 : ',sts.mean(list_moy1_f))
        print('Std (σ) '+type_+'est distance S2 : ',numpy.std(list_moy0_f))
        print('Std (σ) '+type_+'est distance S1 : ',numpy.std(list_moy1_f))
        print('KS test between S1/S2 : ', ' test : ',kolmogorov_smirnov(list_moy0_f,list_moy1_f)[0], ' p-value : ',kolmogorov_smirnov(list_moy0_f,list_moy1_f)[1])
        print('Mann-Whitney rank test  between S1/S2 : ', ' test : ',stats.mannwhitneyu(list_moy0_f,list_moy1_f)[0], ' p-value : ',stats.mannwhitneyu(list_moy0_f,list_moy1_f)[1])

    else : 

        for patient_id in dict_clusters[feature_name][0]:
            dict_keys_f5 = list(all_features[patient_id][feature_name].keys())
            dict_keys_f5.sort()
            if dict_keys_f5 != []:
                list_moy0_f.append(dict_keys_f5[0])

        list_moy0_f = [k*15 for k in list_moy0_f]

        for patient_id in dict_clusters[feature_name][1]:
            dict_keys_f5 = list(all_features[patient_id][feature_name].keys())
            dict_keys_f5.sort()
            if dict_keys_f5 != []:
                list_moy1_f.append(dict_keys_f5[0])

        list_moy1_f = [k*15 for k in list_moy1_f]
        print('Mean f5 smallest distance S2 : ',sts.mean(list_moy0_f))
        print('Mean f5 smallest distance S1 : ',sts.mean(list_moy1_f))
        print('Std (σ) f5 smallest distance S2 : ',numpy.std(list_moy0_f))
        print('Std (σ) f5 smallest distance S1 : ',numpy.std(list_moy1_f))
        print('KS test between S1/S2 : ', ' test : ',kolmogorov_smirnov(list_moy0_f,list_moy1_f)[0], ' p-value : ',kolmogorov_smirnov(list_moy0_f,list_moy1_f)[1])
        print('Mann-Whitney rank test  between S1/S2 : ', ' test : ',stats.mannwhitneyu(list_moy0_f,list_moy1_f)[0], ' p-value : ',stats.mannwhitneyu(list_moy0_f,list_moy1_f)[1])



    
#############################################################################################################################
##################################################### Feature functions #####################################################
#############################################################################################################################

"""
    Finds the distances from lymphocytes to closest cancer cell.
"""

def distances_from_lymphocytes_to_cancer (graph, memoize=True) :
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    
    # Interface cells
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    is_lymphocyte = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"

    # Distances
    distances_lymphocytes_cancer = shortest_distances_between_types(graph, is_lymphocyte, is_cancer)
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, distances_lymphocytes_cancer)
    return distances_lymphocytes_cancer

#############################################################################################################################

"""
    Finds the adjusted distances from lymphocytes to cancer border.
"""

def adjusted_distances_from_lymphocytes_to_cancer_border (graph, memoize=True) :
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    
    # Interface cells
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    is_stroma = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Stroma"
    is_lymphocyte = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"
    interface_triangles = find_interface_triangles(graph, is_cancer, is_stroma)
    border_cancer_cells = set([vertex for vertex in flatten(interface_triangles) if is_cancer(vertex)])
    border_stroma_cells = set([vertex for vertex in flatten(interface_triangles) if is_stroma(vertex)])
    is_border_cancer_cell = lambda vertex : vertex in border_cancer_cells
    is_border_stroma_cell = lambda vertex : vertex in border_stroma_cells

    # Distances
    distances_lymphocytes_cancer = shortest_distances_between_types(graph, is_lymphocyte, is_border_cancer_cell)
    distances_lymphocytes_stroma = shortest_distances_between_types(graph, is_lymphocyte, is_border_stroma_cell)
    distances_adjusted = (numpy.sign([distances_lymphocytes_cancer[i] - distances_lymphocytes_stroma[i] for i in range(len(distances_lymphocytes_cancer))]) * distances_lymphocytes_cancer).tolist()
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, distances_adjusted)
    return distances_adjusted

#############################################################################################################################

"""
    Computes the areas of all triangles at cancer/stroma interface.
"""

def areas_of_triangles_at_interface_cancer_stroma (graph, memoize=True) :
  
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    
    # Interface cells
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    is_stroma = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Stroma"
    interface_triangles = find_interface_triangles(graph, is_cancer, is_stroma)

    # Areas of triangles
    areas = []
    for cluster in interface_triangles :
        for triangle in cluster :
            area = shape_area(triangle, graph.data_frame)
            areas.append(area)
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, areas)
    return areas

#############################################################################################################################

"""
    Computes the density of lymphocyte per band at cancer/stroma interface for a given band size
    Band sizes are fixed at 20µm
"""

def density_lymphocytes (graph,band_size=20, memoize=True) :
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    
    # Interface cells
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    is_stroma = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Stroma"
    is_lymphocyte = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"
    interface_triangles = find_interface_triangles(graph, is_cancer, is_stroma)
    border_cancer_cells = set([vertex for vertex in flatten(interface_triangles) if is_cancer(vertex)])
    border_stroma_cells = set([vertex for vertex in flatten(interface_triangles) if is_stroma(vertex)])
    is_border_cancer_cell = lambda vertex : vertex in border_cancer_cells
    is_border_stroma_cell = lambda vertex : vertex in border_stroma_cells

    # Distances
    distances_stroma_cancer = shortest_distances_between_types(graph, is_stroma, is_border_cancer_cell)
    distances_tumor_cancer = shortest_distances_between_types(graph, is_cancer, is_border_cancer_cell)
    distances_lymphocytes_cancer = shortest_distances_between_types(graph, is_lymphocyte, is_border_cancer_cell)
    distances_lymphocytes_stroma = shortest_distances_between_types(graph, is_lymphocyte, is_border_stroma_cell)
    distances_adjusted = (numpy.sign([distances_lymphocytes_cancer[i] - distances_lymphocytes_stroma[i] for i in range(len(distances_lymphocytes_cancer))]) * distances_lymphocytes_cancer).tolist()
    if not distances_adjusted:
        dict_final={0:[0]}
        if memoize :
            memoize_save(contents_or_target, dict_final)
        return dict_final
    else : 

        distances_tumor_cancer= Convert(distances_tumor_cancer)
        distances_tumor_cancer = [i for i in distances_tumor_cancer if i != 0]
        distance_total = distances_stroma_cancer+distances_tumor_cancer+distances_adjusted



        list_dens_lympho=[]
        list_dens=[]
        dict_bins={}
        dict_bins_lympho={}
        key_to_remove=[]
        dict_final={}

        for dist in distances_adjusted:
            list_dens_lympho.append(int(dist//band_size))
        for dist in distance_total:
            list_dens.append(int(dist//band_size))

        for i in numpy.arange(min(list_dens_lympho),max(list_dens_lympho)+1,1):
            dict_bins_lympho[i]=0
        for l in list_dens_lympho:
            dict_bins_lympho[l]+=1


        for i in numpy.arange(min(list_dens),max(list_dens)+1,1):
            dict_bins[i]=0
        for l in list_dens:
            dict_bins[l]+=1


        for x in dict_bins.keys():
            if not x in dict_bins_lympho.keys():
                key_to_remove.append(x)
        for k in key_to_remove:
            del dict_bins[k]

        for i in dict_bins_lympho.keys():
            if dict_bins_lympho[i]==0:
                dict_final[i]=[0]
            else:
                dict_final[i]=[dict_bins_lympho[i]/dict_bins[i]]
        
        dict_final = {k:v for k,v in dict_final.items() if v != [0]}

        # Memoize and return
        if memoize :
            memoize_save(contents_or_target, dict_final)
        return dict_final
    
#############################################################################################################################

"""
    Computes the ratio of lymphocytes
    Not used
"""

def number_of_lymphocytes(cells):
    number_of_cancer_cells=cells[cells['Category']==3].count()[0]
    number_of_stroma_cells=cells[cells['Category']==2].count()[0]
    number_of_lymphocyte_cells=cells[cells['Category']==1].count()[0]
    number_of_cells=cells.count()[0]
    return number_of_lymphocyte_cells/number_of_cells

#############################################################################################################################

"""
    Computes the number of cancer blobs
    Not used
"""
def number_of_cancer_blobs(graph):
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    return len(find_cell_type_groups(graph,is_cancer))



#############################################################################################################################

"""
    Computes the number of lymphocytes in the biggest blob. Not defined yet
"""

def Feature_Nombre_lymphocytes_dans_tumeur(graph): # NOMBRE DE LYMPHO DANS LA GROSSE TUMEUR
    Hist=[]
    return(Hist)
            
    



#############################################################################################################################

"""
    [0] Computes the number of cells per blobs
    [1] Computes the ratio Number of cells in the biggest cancer blobs / Number of cancer cells.
"""

def number_of_cells_in_cancer_blobs_and_ratio_tumor_new(graph,memoize=True):
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
        
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    blobs=(find_cell_type_groups(graph,is_cancer))
    len_blobs=[]
    for blob in blobs:
        len_blobs.append(len(blob))
    len_blobs.sort()
    ratio_big_tumor=len_blobs[-1]/graph.data_frame.count()[0]
    len_blobs=[ k for k in len_blobs[:-1] if k>0] #drop big blob which is tumort and exclude size1 blob bcz potential error
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, [len_blobs,ratio_big_tumor])
        
    return [len_blobs,ratio_big_tumor]

#############################################################################################################################

"""
    Computes the size of lymphocytes blobs
"""

def number_of_cells_in_lympho_blobs(graph,memoize=True):
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    is_lympho = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"
    blobs=(find_cell_type_groups(graph,is_lympho))
    len_blobs=[]
    for blob in blobs:
        len_blobs.append(len(blob))
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, len_blobs)
    return len_blobs


    
#############################################################################################################################

"""
    Computes the areas of all triangles at cancer/lymphocyte interface.
"""

def areas_of_triangles_at_interface_cancer_lymphocyte(graph, memoize=True) :
  
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target
    
    # Interface cells
    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    is_lymphocyte = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"
    interface_triangles = find_interface_triangles(graph, is_cancer, is_lymphocyte)

    # Areas of triangles
    areas = []
    for cluster in interface_triangles :
        for triangle in cluster :
            area = shape_area(triangle, graph.data_frame)
            areas.append(area)
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, areas)
    return areas
    
    
#############################################################################################################################

"""
    Computes the distance between cancer blobs
"""

def distance_between_cancer_blobs(graph, memoize=True):
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target

    is_cancer = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Cancer"
    blobs=(find_cell_type_groups(graph,is_cancer))

    type_1_vertices = [vertex for vertex in range(graph.vcount()) if is_cancer(vertex)]
    all_distances = all_pairwise_euclidean_distances(graph.data_frame,type_1_vertices,type_1_vertices)

    filtered_df = graph.data_frame[graph.data_frame['Category']==3].reset_index()

    blobs2=[]
    for blob in blobs :
        blob2=[]
        for cell in blob :
            real_index = filtered_df[filtered_df["index"]==cell]
            blob2.append(real_index.index[0])
        blobs2.append(blob2)

    for blob in blobs2 :
        all_distances[blob,blob]=numpy.inf

    distance_final=[]   

    for blob in blobs2 :
        distance_final.append(all_distances[blob,:].min())
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, distance_final)
    return distance_final


#############################################################################################################################

"""
    Computes the distance between lymphocytes blobs
"""

def distance_between_lymphocytes_blobs(graph, memoize=True):
    
    # Result already computed?
    if memoize :
        loaded, contents_or_target = memoize_load(graph)
        if loaded :
            return contents_or_target

    is_lympho = lambda vertex : CATEGORIES[graph.data_frame.iloc[vertex]["Category"]] == "Lymphocyte"
    blobs=(find_cell_type_groups(graph,is_lympho))

    type_1_vertices = [vertex for vertex in range(graph.vcount()) if is_lympho(vertex)]
    all_distances = all_pairwise_euclidean_distances(graph.data_frame,type_1_vertices,type_1_vertices)

    filtered_df = graph.data_frame[graph.data_frame['Category']==1].reset_index()

    blobs2=[]
    for blob in blobs :
        blob2=[]
        for cell in blob :
            real_index = filtered_df[filtered_df["index"]==cell]
            blob2.append(real_index.index[0])
        blobs2.append(blob2)

    for blob in blobs2 :
        all_distances[blob,blob]=numpy.inf

    distance_final=[]   

    for blob in blobs2 :
        distance_final.append(all_distances[blob,:].min())
    
    # Memoize and return
    if memoize :
        memoize_save(contents_or_target, distance_final)
    return distance_final


