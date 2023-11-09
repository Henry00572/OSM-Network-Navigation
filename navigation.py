'''
Navigation through cycle networks (node and basic networks).
'''
__author__ = "Hendrik (Henry572)"
__copyright__ = "No rights reserved"
__credits__ = ["Hendrik (Henry572)"]
__license__ = "CC0 1.0"
__version__ = "1.0"

import json
import geopy.distance

# General utility functions

def calculate_distance(lon1, lat1, lon2, lat2):
    '''Calculates the distance between the given coordinates

    Arguments:
    lon1:   Longitude of the first point (E/W)
    lat1:   Latitude of the first point (N/S)
    lon2:   Longitude of the second point (E/W)
    lat2:   Latitude of the second point (N/S)
    Return values:
    The distance between the points in meters
    Side effects:
    None
    Dependencies:
    geopy

    The calculation using Geopy considers a geodesic on the WGS-84 ellipsoid
    '''
    return geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).meters

def restructure_list(object_list, object_type=None, keep_type=False, keep_order=False):
    '''Restructures a given list of objects into a dictionary using the object "id" as the key
    
    Arguments:
    object_list:    A list of dictionaries each having at least the key "id"
    object_type:    objects will be filtered so only those having a "type"=object_type tag are returned
                    Optional, if none, all objects are returned
    keep_type:      Whether to keep the key "type" in the dictionary
                    Optional, default ist False
    keep_position:  Whether to insert a tag "position" into each object
                    that describes the position in the object_list
    Return values:
    A dictionary where the keys are the "id" of the objects wich are themselves dictionaries
    Side effects:
    None
    '''
    object_dict = {}
    i = 0
    for obj in object_list:
        if not(object_type) or obj["type"] == object_type:
            object_dict[obj["id"]] = {key: obj[key] for key in obj
                                      if key != "id" and (keep_type or key != "type")}
            if keep_order:
                object_dict[obj["position"]] = i
        i += 1
    return object_dict

def combine_directions(directions):
    '''Combines the possible directions that a way can be traversed along
    
    Arguments:
    directions:     A list of directions or roles, i.e. "", "forward" or "backward", duplicates possible
    Return values:
    A string representing the combined effect of the possible directions
    If "" or both "forward" and "backward" are given, then the way can be used in both directions
    Otherwise only in one direction
    Side effects:
    None
    '''
    if "" in directions or ("forward" in directions and "backward" in directions):
        return ""
    else:
        return directions[0]

def add_route_data(route, coordinate_dict, add_distances=True):
    '''Adds location data to a route given by a list of node ids

    Arguments:
    route:              A list of node ids describing the route
    coordinate_dict:    A dictionary mapping node ids to dictionaries with "lat" and "lon" data
    add_distances:      Whether to add a cumulative distance from the start to each node
                        Optional, the default ist True
    Return values:
    A list of dictionaries each having at least an "id", "lat" and "lon" and possibly "distance" tag
    The order of the route is preserved
    Side effects:
    None
    '''
    for i in range(len(route)):
        route[i] = {"id": route[i]} if add_distances else route[i]
        route[i]["lon"] = coordinate_dict[route[i]["id"]]["lon"]
        route[i]["lat"] = coordinate_dict[route[i]["id"]]["lat"]
    if add_distances:
        for i in range(len(route)):
            route[i]["distance"] = (route[i-1]["distance"]
                                    +calculate_distance(coordinate_dict[route[i-1]["id"]]["lon"], coordinate_dict[route[i-1]["id"]]["lat"],
                                                        coordinate_dict[route[i]["id"]]["lon"], coordinate_dict[route[i]["id"]]["lat"])
                                    if i != 0 else 0)
    return route

## Using topology

def get_child_nodes_of_way(way_id, way_dict):
    '''Returns the nodes of the way
    
    Arguments:
    way_id:         The id of the way
    way_dict:       A dictionary of ways with members
                    It only has to contain the specified way
    Return values:
    A list of node ids in the order they have in the way
    Side effects:
    None
    '''
    return way_dict[way_id]["nodes"]

def get_parent_ways_of_node(node_id, relation_id, way_dict, relation_dict):
    '''Returns the ways in the route relation that contain the given node
    
    Arguments:
    node_id:            The id of the node
    relation_id:        The id of the relation
    way_dict:           A dictionary of ways with members
                        It only has to contain the member ways of the relation
    relation_dict:      A dictionary of relations with members
                        It only has to contain the relation with given relation_id
    Return values:
    A list of way ids corresponding to the ways that contain the node
    The order of the ways will be the same as their order in the relation
    Side effects:
    None
    '''
    return [member["ref"] for member in relation_dict[relation_id]["members"]
            if member["type"] == "way" and node_id in way_dict[member["ref"]]["nodes"]]

def get_neighbours_of_node_in_way(node_id, way_id, way_dict):
    '''Returns a list of the predeccessor and successor of the node in the way, if they exist
    
    Arguments:
    node_id:        The id of the node
    way_id:         The id of the way that has to contain the node
    way_dict:       A dictionary of ways with members
                    It only has to contain the specified way
    Return values:
    A list of either one or two neighbour node ids in the order they appear in the way
    Side effects:
    None
    '''
    way_nodes = get_child_nodes_of_way(way_id, way_dict)
    position = way_nodes.index(node_id)
    if position == 0:
        return [way_nodes[1]]
    elif position == len(way_nodes)-1:
        return [way_nodes[position-1]]
    else:
        return [way_nodes[position-1], way_nodes[position+1]]

def get_members_of_relation(relation_id, member_type, relation_dict):
    '''Returns the members of the relation with given type
    
    Arguments:
    relation_id:        The relations's id
    member_type:        The type of the members that should be returned
    relation_dict:      A dictionary of relations with members
                        It only has to contain the relation with given relation_id
    Return values:
    A list of ids referring to the selected members in the order they have in the relation
    Side effects:
    None
    '''
    return [member["ref"] for member in relation_dict[relation_id]["members"] if member["type"] == member_type]

def get_roles_of_member_in_relation(relation_id, member_id, member_type, relation_dict):
    '''Returns all the roles that the given member has in the relation
    
    Arguments:
    relation_id:        The id of the relation
    member_id:          The id of the member
    member_type:        The type of the member
    relation_dict:      A dictionary of relations with members
                        It only has to contain the relation with given relation_id
    Return values:
    A list of all the roles the given member has in the relation
    in the order they appear in the relation, duplicates are possible
    Side effects:
    None
    '''
    return [member["role"] for member in relation_dict[relation_id]["members"]
            if member["ref"] == member_id and member["type"] == member_type]

def get_parent_relations_of_intersection(intersection_id, intersection_dict):
    '''Returns all relations containing the given node in one of their member ways
    
    Arguments:
    node_id:                The id of the node
    intersection_dict:      A dictionary of intersections with parent relations
                            It only has to contain the intersection with given intersection_id
    Return values:
    A list of relation ids that have ways as members containing the node
    The ids will be in the same order as in the intersection_dict
    Side effects:
    None
    '''
    try:
        return intersection_dict[intersection_id]["routes"]
    except KeyError as e:
        print("Fehler:", e)
        print(len(intersection_dict.keys()))

def get_parent_relations_of_intersections(intersection_dict, intersection_ids=None):
    '''Returns all relations containing one of the given nodes in one of their member ways
    
    Arguments:
    intersection_dict:      A dictionary of intersections with parent relations
                            It only has to contain the intersection with given intersection_id
    intersection_ids:       A list of ids of intersections
                            Optional, if None, all intersections in intersection_dict will be considered
    Return values:
    A list of relation ids that have ways as members containing intersections in the list
    Side effects:
    None
    '''
    if not(intersection_ids):
        intersection_ids = list(intersection_dict.keys())
    return list(set().union(*(get_parent_relations_of_intersection(intersection_id, intersection_dict) for  intersection_id in intersection_ids)))

## IO

def load_data(topology_filename, geometry_filename):
    '''Loads the neccessary data from files
    
    Arguments:
    topology_filename:      The file containing routes and ways with members
                            including its path and/or file extension
    geometry_filename:      The file contaning nodes with coordinates
                            including its path and/or file extension
    Return values:
    topology, geometry:     All list of dictionaries containing the OSM objects
    Side effects:
    None
    Dependencies:
    json

    The files with appropiate format can be created by using the overpass queries
    and exporting the result as raw OSM data
    '''
    f = open(topology_filename)
    topology_load = json.load(f)["elements"]
    f.close()
    f = open(geometry_filename)
    geometry_load = json.load(f)["elements"]
    f.close()
    return topology_load, geometry_load

def export_route(route, file_name, file_format, override=True):
    '''Exports the route to a file
    
    Arguments:
    route:              A list of dictionaries with node ids and location data
                        Currently, "id", "lon" and "lat" are required
                        and "distance" is optional
    file_name:          The name or path for the target file without a file extension
    file_format:        The format used for the export
                        Currently supported is "csv"
    override:           Whether to override the file, otherwise append to it
                        Optional, the default is True
    Return value:
    None
    Side effects:
    Writes the coordinates of the given points to the specified file and may override any existing file content
    For "csv": Each node will have one line with "id", "lon", "lat" and possibly "distance" separated by ", "
    May print out an error if a different file format is given or throw an error if data is missing
    Dependencies:
    json
    '''
    write_mode = {True: "w", False: "a"}
    match file_format:
        case "csv":
            f = open(str(file_name) + ".csv", write_mode[override])
            for node in route:
                f.write(str(node["id"]) +
                        "," + str(node["lon"]) +
                        "," + str(node["lat"]) +
                        ("," + str(node["distance"]) if "distance" in node else "") +
                        "\n")
            f.close()
        case "json":
            f = open(str(file_name) + ".json", write_mode[override])
            json.dump(route, f, indent=4)
            f.close()
        case other:
            print("Invalid file format")

def import_network_adjacency(file_name):
    '''Imports a network adjacency dictionary from a json file
    
    Arguments:
    filename:               The path/name of the file, ideally in ".json" format
    Return values:
    An adjacency dictionary for the network intersections
    Side effects:
    None

    Node ids that are used as dictionary keys have to be converted from strings to ints
    because json only allows for strings as keys, but this script uses ints
    '''
    f = open(file_name)
    network_adjacency = json.load(f)
    f.close()
    return {int(intersection_id): {int(neighbour_id): distance
                                   for neighbour_id,distance in intersection_adjacency.items()}
            for intersection_id, intersection_adjacency in network_adjacency.items()}

def export_network_adjacency(network_adjacency, file_name):
    '''Exports the calculated network adjacency to a file
    
    Arguments:
    network_adjacency:      Adjacency dictionary as returned by network_adjacency()
    filename:               The path/name of the file, ideally in ".json" format
    Return values:
    None
    Side effects:
    Writes to the specified file, will override existing content
    Dependencies:
    json
    '''
    f = open(file_name, "w")
    json.dump(network_adjacency, f, indent=4)
    f.close()

# Specific algorithms

def generate_intersection_list(way_dict, relation_dict):
    '''Generates the intersection dictionary with nodes where different routes join
    
    Arguments:
    relation_dict:          A dictionary of relations and their members
    way_dict:               A dictionary of ways and their members
    Return values:
    A dictionary with all the intersections' ids as keys and
    a list of ids of relations leading through the intersection as values
    Side effects:
    None

    Intersections are nodes where a route might change between relations
    They're at the start/end of a relation
    or where two relations take different ways to/from the node
    Relations running parallel for some ways will only lead to intersections
    at the start and end of the parallel section
    '''
    all_nodes = {}
    intersection_dict = {}
    # Recurse down through all relations and add the parent relations to each node
    for relation_id in relation_dict:
        for way_id in get_members_of_relation(relation_id, "way", relation_dict):
            for node_id in get_child_nodes_of_way(way_id, way_dict):
                if node_id in all_nodes:
                    all_nodes[node_id].add(relation_id)
                else:
                    all_nodes[node_id] = {relation_id}
    all_nodes = {node: list(routes) for node,routes in all_nodes.items()}
    # Iterate over the nodes and check if they are an intersection
    for node_id in all_nodes:
        parent_ways_per_relation = [get_parent_ways_of_node(node_id, relation_id, way_dict, relation_dict)
                                    for relation_id in all_nodes[node_id]]
        parent_ways = list(set([way for ways in parent_ways_per_relation for way in ways]))
        neighbour_nodes = list(set([neighbour for way_id in parent_ways
                                    for neighbour in get_neighbours_of_node_in_way(node_id,way_id,way_dict)]))
        # Core logic for identifying intersections
        if (len(neighbour_nodes) >= 3
            or any(set(parent_ways) != set(parent_ways_per_relation[0]) for parent_ways in parent_ways_per_relation)):
            intersection_dict[node_id] = {"routes": all_nodes[node_id]}
    return intersection_dict

def generate_relation_adjacency(relation_id, way_dict, relation_dict, coordinate_dict):
    '''Creates an adjacency dictionary for the route relation with given id
    
    Arguments:
    relation_id:        The id of the route relation
    way_dict:           A dictionary of ways with members
                        It only has to contain the member ways of the relation
    relation_dict:      A dictionary of relations with members
                        It only has to contain the relation with given relation_id
    coordinate_dict:    A dictionary of dictionaries where the node ids can be used as keys
                        and "lat" and "lon" as subkeys for latitute and longitude of a point
    Return values:
    A dictionary representing the adjacency of nodes in the relation's member ways
    The keys are node ids and their values are dictionaries with node ids as keys and
    a dictionary of distance and nodes as value
    Side effects:
    None

    All nodes of ways in the relation are considered. For each node,
    its parent ways with their roles describing the possible directions are considered and the next nodes
    in those directions are saved together with the distance.
    For a linear relation this creates a doubly linked list, while the result for a relation with forward/backward ways
    will lead to more successors at the branching nodes
    '''
    # Find all the nodes of ways in the relation
    nodes = []
    for way_id in get_members_of_relation(relation_id, "way", relation_dict):
        nodes += get_child_nodes_of_way(way_id, way_dict)
    nodes = list(set(nodes))
    relation_adjacency = {}
    for node in nodes:
        relation_adjacency[node] = {}
        parent_ways = get_parent_ways_of_node(node, relation_id, way_dict, relation_dict)
        for parent_way in parent_ways:
            way_nodes = get_child_nodes_of_way(parent_way, way_dict)
            index = way_nodes.index(node)
            direction = combine_directions(get_roles_of_member_in_relation(relation_id, parent_way, "way", relation_dict))
            # Check in which direction the parent way can be used and calculate the distance to the successor in that direction
            if direction in ["", "forward"] and index != len(way_nodes)-1:
                next_node = way_nodes[index+1]
                relation_adjacency[node][next_node] = {"distance": calculate_distance(coordinate_dict[node]["lon"],
                                                                                      coordinate_dict[node]["lat"],
                                                                                      coordinate_dict[next_node]["lon"],
                                                                                      coordinate_dict[next_node]["lat"]),
                                                       "nodes": [node, next_node]}
            if direction in ["", "backward"] and index != 0:
                next_node = way_nodes[index-1]
                relation_adjacency[node][next_node] = {"distance": calculate_distance(coordinate_dict[node]["lon"],
                                                                                      coordinate_dict[node]["lat"],
                                                                                      coordinate_dict[next_node]["lon"],
                                                                                      coordinate_dict[next_node]["lat"]),
                                                       "nodes": [node, next_node]}
    return relation_adjacency

def generate_network_adjacency(intersection_dict, way_dict, relation_dict, coordinate_dict):
    '''Creates an adjacency dictionary where nodes connected by a relation are adjacent
    
    Arguments:
    intersection_dict:  A list of node ids that form the intersections of the network
    way_dict:           A dictionary of ways with members
                        It only has to contain the member ways of the relation
    relation_dict:      A dictionary of relations with members
                        It only has to contain the relation with given relation_id
    coordinate_dict:    A dictionary of dictionaries where the node ids can be used as keys
                        and "lat" and "lon" as subkeys for latitute and longitude of a point
    Return values:
    A dictionary representing the adjacency of nodes in the relation's member ways
    The keys are node ids and their values are dictionaries with node ids as keys and
    a dictionary containing "distance" and "route"
    Side effects:
    Prints some status updates to standard output

    Adjacent nodes will have a distance according to their minimum distance along their shared relations.
    This may differ between the two directions A->B and B->A.
    The "route" of two nodes is a list of node ids describing the shortest path, including start and end.
    '''
    print("Creating network adjacency")
    network_adjacency = {}
    parent_relations = get_parent_relations_of_intersections(intersection_dict, list(intersection_dict.keys()))
    relation_adjacencies = {}
    # Track progress
    current_i = 0
    total_i = len(parent_relations)
    # For each relation, create its adjacency and store them
    for relation_id in parent_relations:
        print("Creating relation adjacency for relation", relation_id, "(", current_i, "/", total_i, ")")
        current_i += 1
        relation_adjacencies[relation_id] = generate_relation_adjacency(relation_id, way_dict, relation_dict, coordinate_dict)
    print("Combining relation adjacencies")
    # Track progress
    current_i = 0
    total_i = len(intersection_dict)
    # Consider each pair of intersections
    for intersection_id in intersection_dict:
        print("Considering intersection", intersection_id, "(", current_i, "/", total_i, ")")
        current_i += 1
        network_adjacency[intersection_id] = {}
        for other_intersection_id in intersection_dict:
            common_parent_relations = list(set(get_parent_relations_of_intersection(intersection_id, intersection_dict))
                                         & set(get_parent_relations_of_intersection(other_intersection_id, intersection_dict)))
            if common_parent_relations:
                route_options = [dijkstra(intersection_id, other_intersection_id, relation_adjacencies[common_parent_relation])
                                 for common_parent_relation in common_parent_relations]
                route_options = [{"distance": route_option[0],
                                  "nodes": fill_intermediate_nodes_in_route(other_intersection_id,
                                                                            route_option[1],
                                                                            relation_adjacencies[common_parent_relation])}
                                 for route_option,common_parent_relation in zip(route_options,common_parent_relations)]
                shortest_distance = min([route_option["distance"] for route_option in route_options])
                if shortest_distance != float("inf"):
                    shortest_route = [route_option["nodes"] for route_option in route_options
                                      if route_option["distance"] == shortest_distance][0]
                    network_adjacency[intersection_id][other_intersection_id] = {"distance": shortest_distance,
                                                                                 "nodes": shortest_route}
                else:
                    network_adjacency[intersection_id][other_intersection_id] = {"distance": float("inf")}
    return network_adjacency

def dijkstra(start_id, end_id, adjacency):
    '''Calculates the shortest path between the start and end nodes using Dijkstra's algorithm
    
    Arguments:
    start_id:       id of the start node
    end_id:         id of the end node
    adjacency:      dictionary of dictionaries where the keys are node ids
                    and the inner dictionaries contain distances and possibly node lists for the distance
    Return values:
    The calculated distance from start to end
    and a dictionary mapping each node id out of adjacency to a dictionary with
    cumulative "distance" from the start and each node's "predeccessor"
    Side effects:
    None

    The algorithm makes use of the hierarchy of route relations forming a network.
    If adjacency is the output of generate_relation_adjacency, all nodes will be considered and arbitrary
    routing along the relation is possible.
    If adjacency is the output of generate_network_adjacency, only intersections will be considered.
    Their distance to other intersections in the same route has already been computed
    so all intermediate nodes can now be ignored.
    '''
    distances = {start_id: {"distance": 0, "predeccessor": None}}
    if start_id == end_id:
        return 0, distances
    current = start_id
    unvisited = list(adjacency.keys())
    # Run as long as some nodes/intersections are unvisited
    while len(unvisited) > 0:
        for neighbour in adjacency[current]:
            # Update the distances of the current's neighbours
            if ((neighbour in distances and distances[neighbour]["distance"]
                                           > distances[current]["distance"]
                                           + adjacency[current][neighbour]["distance"])
                or neighbour not in distances):
                distances[neighbour] = {"distance": distances[current]["distance"] + adjacency[current][neighbour]["distance"],
                                        "predeccessor": current}
        unvisited.remove(current)
        # possible_next contains the unvisited node(s) closest to the start or nothing
        possible_next = [node for node in list(set(unvisited) & set(distances.keys()))
                         if distances[node]["distance"] == min([distances[n]["distance"]
                                                                for n in list(set(unvisited) & set(distances.keys()))
                                                                if distances[n]["distance"] > 0])]
        if len(possible_next) > 0:
            current = possible_next[0]
        else:
            break
    if end_id not in distances:
        distances[end_id] = {"distance": float("inf")}
    return distances[end_id]["distance"], distances

def fill_intermediate_nodes_in_route(end_id, distances, adjacency):
    '''Creates a route including non-intersection nodes

    Argumemts:
    distances:      The distances calculated by dijkstra
                    A dictionary mapping intersection ids to a dictionary
                    of "distance" from the start and "predeccessor" towards the start
    adjacency:      A dictionary of dictionaries where the keys are node ids
                    and the inner dictionaries contain distances and node lists for the distance
    Return values:
    A list of all node's ids the route passes through
    Side effects:
    None

    The distances form a singly linked list from end_id to the start which can easily be followed.
    While doing this, the intermediate nodes between intersections 
    '''
    if len(distances) == 1:
        return list(distances.keys())
    if distances[end_id]["distance"] == float("inf"):
        return None
    route_nodes = [end_id]
    current_last = end_id
    while distances[current_last]["distance"] != 0:
        route_nodes = ([distances[current_last]["predeccessor"]]
                       + adjacency[distances[current_last]["predeccessor"]][current_last]["nodes"][1:-2]
                       + route_nodes)
        current_last = route_nodes[0]
    return route_nodes


TOPOLOGY_FILENAME = "data/all_data.json"
GEOMETRY_FILENAME = "data/all_points.json"

TOPOLOGY, GEOMETRY = load_data(TOPOLOGY_FILENAME, GEOMETRY_FILENAME)

ALL_WAYS = restructure_list(TOPOLOGY, object_type="way", keep_type=True)
ALL_RELATIONS = restructure_list(TOPOLOGY, object_type="relation")
ALL_COORDINATES = restructure_list(GEOMETRY, object_type="node")
ALL_INTERSECTIONS = generate_intersection_list(ALL_WAYS, ALL_RELATIONS)

#start = 313084340 # KP 55
#end = 2367715463 # S Eller Süd S
#end = 246701723  # Unterbilk
#start = 715750570 # KP 63
#end = 1567989038 # Südlich KP 89
#start = 274353880 # Östlich KP 89
#start = 6988366473 # Karolingerstraße
#end = 2719301824 # KP 89
end = 3446341126 # Rothenbergstr/Kleiner Torfbruch
#end =  257766855 # Waldspielplatz
#end = 11185300932 # Im Broich
#end = 258087981 # KP 05
#end = 3585257573 # Am Schönenkamp/Altenbrückstr
#end = 1544063213 # Hassels Kirche
#start = 257555533 # Abzweig Botanischer Garten S
#end = 241783117 # Lessingplatz
#start = 248908184 # Hellerhof, Mühlenbach
#start = 60601938 # Urdenbacher Allee/Koblenzer Str
#start = 38417819 # Urdenbacher Allee/Haydnstr
start = 256329617 # Bonner Str/Hügelstr
#end = 5060973542 # Werstener Feld/Opladener Str
#end = 1615052289 # Werstener Feld/Ickerswarder Str SO
#end = 253810569 # Kö
#end = 241429117 # KP 71

network_adjacency = generate_network_adjacency(ALL_INTERSECTIONS, ALL_WAYS, ALL_RELATIONS, ALL_COORDINATES)
export_network_adjacency(network_adjacency, "output/network_adjacency.json")
#network_adjacency = import_network_adjacency("network_adjacency.json")
print("Starting dijkstra on network")
distance, distances = dijkstra(start, end, network_adjacency)
route = fill_intermediate_nodes_in_route(end, distances, network_adjacency)
print("Navigation result:", distance)
print(route)
export_route(add_route_data(route, ALL_COORDINATES), "output/route_coordinates_" + str(start) + "_" + str(end), "csv")

#export_route(add_route_data(intersection_list, location_list), "output/intersections", "csv")
