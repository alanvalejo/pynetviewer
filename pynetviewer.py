#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PyNetViewer: A tool for visualization of bipartite, k-partite and heterogeneous networks

Copyright (C) 2017 Alan Valejo <alanvalejo@gmail.com> All rights reserved

PyNetViewer is a python and igraph based tool for visualization of bipartite, k-partite and heterogeneous networks.
The main aim of the PyNetViewer is the visualization of the benchmark networks synthesized by the Bnoc tool.

This program comes with ABSOLUTELY NO WARRANTY. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS
WITH YOU.

Owner or contributors are not liable for any direct, indirect, incidental, special, exemplary, or consequential
damages, (such as loss of data or profits, and others) arising in any way out of the use of this software,
even if advised of the possibility of such damage.

This program is free software and distributed in the hope that it will be useful: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with this program. If not,
see http://www.gnu.org/licenses/.

Giving credit to the author by citing the papers [1]

[1] Valejo, Alan and Goes, F. and Romanetto, L. M. and Oliveira, Maria C. F. and Lopes, A. A., A benchmarking tool
for the generation of bipartite network models with overlapping communities, in Knowledge and information systems,
p. 1-29, 2019 doi: https://doi.org/10.1007/s10115-019-01411-9

Warning: The original implementation (i.e. paper version [1]) is deprecated. This software is a new version, more robust
and fast. There may be divergences between this version and the original algorithm. If you looking for the original
version used in the paper don't hesitate to contact the authors.
"""
import sys

import igraph
import os
import numpy
import random
import inspect
import math

from colour import Color
from fa2 import ForceAtlas2

import models.args as args
import models.helper as helper
import models.helperigraph as helperigraph

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1.0'
__date__ = '2020-04-25'

def build_layout(graph, options):
    coords = []
    w = 0.6 * min(numpy.asarray(options.bbox))
    r = 0.8 * w / 2
    n = len(options.vertices)

    teta = numpy.linspace(-math.pi / 2, 2 * math.pi - math.pi / 2, num=n + 1)
    p = [numpy.asarray([r * math.cos(t), r * math.sin(t)]) for t in teta]

    for layer in range(n):
        p0 = p[layer]
        p1 = p[layer + 1]
        ni = options.vertices[layer]
        vec = p1 - p0
        s = 0.4 / (ni - 1)
        ci = [p0 + (0.3 + s * i) * vec for i in range(ni)]
        coords += ci

    layout = igraph.Layout(coords)

    return layout

def giant(graph):
    components = graph.components()
    components_sizes = components.sizes()
    giant_component_index = components_sizes.index(max(components_sizes))
    return components.subgraph(giant_component_index)

def number_of_components(graph):
    components = graph.components()
    components_sizes = components.sizes()
    return len(components_sizes)

def compute_layout(graph, options, boundary_edges=None):
    if boundary_edges:
        gcopy = graph.copy()
        gcopy.delete_edges(boundary_edges)
    else:
        gcopy = graph
    if igraph.__version__ != '0.8.0':
        if options.layout_name == 'heterogeneous':
            graph.vs['layout'] = build_layout(graph, options)
            for edge in graph.es():
                u, v = edge.tuple
                type_u, type_v = graph.vs['type'][u], graph.vs['type'][v]
                edge['edge_curved'] = 0.2 if type_u != type_v else -2.0
        else:
            graph.vs['layout'] = gcopy.layout(options.layout_name)
    else:

        if options.layout_name == 'forceatlas2':
            forceatlas2 = ForceAtlas2(
                # Behavior alternatives
                outboundAttractionDistribution=options.layout_hub_attraction  # Dissuade hubs
                , linLogMode=False  # NOT IMPLEMENTED
                , adjustSizes=False  # Prevent overlap NOT IMPLEMENTED
                , edgeWeightInfluence=1.0
                # Performance
                , jitterTolerance=10.0
                , barnesHutOptimize=True
                , barnesHutTheta=1.2
                , multiThreaded=False  # NOT IMPLEMENTED
                # Tuning
                , scalingRatio=options.layout_scaling_ratio
                , strongGravityMode=False
                , gravity=options.layout_gravity
                # Log
                , verbose=True
            )
            graph.vs['layout'] = forceatlas2.forceatlas2_igraph_layout(
                gcopy
                , pos=None
                , iterations=options.layout_niter
                , weight_attr="weight"
            )
        if options.layout_name == 'fr':
            graph.vs['layout'] = gcopy.layout(options.layout_name, niter=options.layout_niter)
        elif options.layout_name == 'graphopt':
            graph.vs['layout'] = gcopy.layout(options.layout_name, niter=options.layout_niter)
        elif options.layout_name == 'lgl':
            graph.vs['layout'] = gcopy.layout(options.layout_name, coolexp=1.5, repulserad=-1, cellsize=-1)
        elif options.layout_name == 'heterogeneous':
            graph.vs['layout'] = build_layout(graph, options)
            for edge in graph.es():
                u, v = edge.tuple
                type_u, type_v = graph.vs['type'][u], graph.vs['type'][v]
                edge['edge_curved'] = 0.2 if type_u != type_v else -2.0
        else:  # drl, sugiyama, auto, circle, random, rt, rt_circular, mds, star, dh, bipartite, kk
            graph.vs['layout'] = gcopy.layout(options.layout_name)

        if options.layout_to_radial:
            graph.vs['layout'].to_radial(min_angle=50, max_angle=10, min_radius=0.0, max_radius=1.0)

    with open(options.output + '.layout', 'w+') as f:
        for xy in graph.vs['layout']:
            f.write(str(xy[0]) + ',' + str(xy[1]) + '\n')

def compute_min_max(vector, new_min, new_max):
    result = vector.copy()
    old_min, old_max = min(result), max(result)
    for key, item in enumerate(vector):
        new_value = helper.remap(item, old_min, old_max, new_min, new_max)
        if new_value is None or math.isnan(new_value):
            new_value = new_min
        result[key] = new_value
    return result

def get_boundary_edges(graph):
    boundary_edges = []
    for edge in graph.es():
        if graph.vs[edge.tuple[0]]['vertex_color'] != graph.vs[edge.tuple[1]]['vertex_color']:
            boundary_edges.append(edge)
    return boundary_edges

def read_membership(graph, file_membership):
    graph['comms'] = [0] * graph['layers']
    graph['overlapping'] = []
    if options.file_membership:
        with open(options.file_membership, 'r') as f:
            for vertex, comms in enumerate(f):
                members = set(map(int, comms.split()))
                if len(members) > 1:
                    graph['overlapping'].append(vertex)
                graph.vs[vertex]['membership'] = members
        for layer in range(graph['layers']):
            vertices = graph.vs.select(type=layer)['index']
            graph['comms'][layer] = max(list(set().union(*graph.vs[vertices]['membership']))) + 1

def read_weight(graph, file):
    graph.vs['weight'] = numpy.loadtxt(file)

def read_layout(graph, file):
    graph.vs['layout'] = [None] * graph.vcount()
    with open(file, 'r') as f:
        for vertex, line in enumerate(f):
            item = list(map(float, line.strip().split(',')))
            graph.vs[vertex]['layout'] = [item[0], item[1]]

def community_detection(graph, algorithm="fastgreedy", k=3, output='output_file'):
    k_componentes = number_of_components(graph)
    if k_componentes > k:
        log.warning('The number of components in the graph is less than the number of communities defined.')
        log.warning('Number of components: ' + str(k_componentes))
        log.warning('Number of communities defined: ' + str(k))
        sys.exit(1)
    if algorithm == "fastgreedy":
        cl = graph.community_fastgreedy(weights='weight')
        membership = cl.as_clustering(k).membership
        for vertex, members in enumerate(membership):
            graph.vs[vertex]['membership'] = set([members])
        graph['comms'] = [k] * graph['layers']
        with open(output + '.membership', 'w+') as f:
            for vertex in graph.vs():
                f.write(' '.join(map(str, list(vertex['membership']))) + '\n')
    else:
        parser.error('There are no ' + str(algorithm) + ' algorithm.')

def compute_shapes(graph):
    types = ['rectangle', 'circle', 'triangle-up', 'triangle-down']
    while len(types) < graph['layers']:
        types += types
    for layer in range(graph['layers']):
        for vertex in graph.vs.select(type=layer):
            vertex['vertex_shape'] = types[layer]

def read_color(graph, file):
    with open(file, 'r') as f:
        for index, line in enumerate(f):
            graph.vertex['color'].append(line.strip())

def compute_vertex_color_by_membership(graph, eq_colors=False, output='output_file'):
    colors = []
    if eq_colors:
        total_colors = max(graph['comms'])
    else:
        total_colors = sum(graph['comms'])

    for i in range(0, total_colors):
        colors.append('#' + '%06X' % random.randint(0, 0xFFFFFF))

    with open(output + '.color', 'w+') as f:
        for color in colors:
            f.write('\n'.join(color))

    for layer in range(graph['layers']):
        for vertex in graph.vs.select(type=layer):
            member = next(iter(vertex['membership']))
            if eq_colors:
                vertex['vertex_color'] = colors[member]
            else:
                vertex['vertex_color'] = colors[member + sum(graph['comms'][:layer])]

if __name__ == '__main__':

    # Setup parse options command line
    current_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    parser = args.setup_parser(current_path + '/args/pynetviewer.json')
    options = parser.parse_args()
    args.update_json(options)
    args.check_output(options)

    # Log instanciation
    log = helper.initialize_logger(dir='log', output='log')

    # Check required fields
    if options.input is None:
        parser.error('required -f [input] arg.')
    if options.vertices is None:
        parser.error('required -v [number of vertices for each layer] arg.')

    graph = helperigraph.load(options.input, options.vertices, type_filename=options.file_type)

    graph.vs['membership'] = None
    graph['overlapping'] = None
    graph['comms'] = None
    graph['overlapping'] = None
    graph.vs['weight'] = 1
    graph.vs['vertex_color'] = 'black'
    graph.vs['vertex_size'] = options.vertex_size_min
    graph.vs['vertex_shape'] = graph.vcount() * ['circle']
    graph.vs['vertex_frame_color'] = options.vertex_frame_color
    graph.vs['vertex_frame_width'] = options.vertex_frame_width
    boundary_edges = None

    # Read files
    if options.file_membership:
        read_membership(graph, options.file_membership)
    if options.file_weight:
        read_weight(graph, options.file_weight)
    if options.file_layout:
        read_layout(graph, options.file_layout)
    if options.file_color:
        read_color(graph, options.file_color)

    graphs = [graph]

    if options.split_projections:
        graphs = []
        projections = graph.bipartite_projection(graph.vs['type'])
        for key in options.split_projections:
            graphs.append(projections[key])

    options.original_output = options.output
    for key, graph in enumerate(graphs):

        if options.split_projections:
            options.output = options.original_output + '-p-' + str(key)

        # Filter graph
        if options.only_giant_component:
            graph = giant(graph)
        if options.delete_edge_by_weight_le:
            for weight in options.delete_edge_by_weight_le:
                graph.es.select(weight_le=weight).delete()
        if options.delete_vertex_by_degree_le:
            for degree in options.delete_vertex_by_degree_le:
                graph.vs.select(_degree=degree).delete()

        # Community detection
        if options.community_detection_algorithm:
            community_detection(graph, algorithm=options.community_detection_algorithm,
                                k=options.number_of_communities, output=options.output)

        # Style
        if options.vertex_size_by == 'degree':
            graph.vs['weight'] = compute_min_max(graph.degree(), options.vertex_size_min,
                                                 options.vertex_size_max)
        if options.vertex_size_by == 'strength':
            graph.vs['weight'] = compute_min_max(graph.strength(weights='weight'), options.vertex_size_min,
                                                 options.vertex_size_max)
        if options.vertex_size_by == 'weight':
            graph.vs['weight'] = compute_min_max(graph.vs['weight'], options.vertex_size_min, options.vertex_size_max)

        if options.vertex_color_by == 'weight':
            unique_degrees = list(numpy.unique(graph.strength()))
            blue = Color("#000000")
            red = Color("#CC7B7B")
            colors = list(blue.range_to(red, len(unique_degrees)))
            graph.vs['vertex_color'] = [str(colors[unique_degrees.index(degree)]) for degree in graph.strength()]
        elif options.vertex_color_by == 'degree':
            unique_weights = list(numpy.unique(graph.vs['weight']))
            blue = Color("#000000")
            red = Color("#CC7B7B")
            colors = list(blue.range_to(red, len(unique_weights)))
            graph.vs['vertex_color'] = [str(colors[unique_weights.index(weight)]) for weight in graph.vs['weight']]
        elif options.vertex_color_by == 'strength':
            unique_weights = list(numpy.unique(graph.strength(weights='weight')))
            blue = Color("#000000")
            red = Color("#CC7B7B")
            colors = list(blue.range_to(red, len(unique_weights)))
            graph.vs['vertex_color'] = [str(colors[unique_weights.index(weight)]) for weight in graph.strength(weights='weight')]
        elif options.vertex_color_by == 'membership':
            if not options.community_detection_algorithm and not options.file_membership:
                log.warning('To coloring vertex by its membership provide a membership file ou set a community '
                            'detection algorithm.')
                sys.exit(1)
            compute_vertex_color_by_membership(graph, eq_colors=options.eq_colors, output=options.output)

        if options.overlapping_paint and graph['overlapping']:
            overlapping = graph['overlapping']
            graph.vs[overlapping]['vertex_color'] = str(options.overlapping_color)
            graph.vs[overlapping]['vertex_frame_color'] = str(options.vertex_frame_color)
            graph.vs[overlapping]['vertex_shape'] = str(options.overlapping_shape)

        compute_shapes(graph)
        graph.es['weight'] = compute_min_max(graph.es['weight'], options.edge_weight_min, options.edge_weight_max)
        graph.es['opacity'] = compute_min_max(graph.es['weight'], options.edge_opacity_min, options.edge_opacity_max)
        for edge in graph.es():
            edge['opacity'] = "rgba(0,0,0," + str(edge['opacity']) + ")"

        if options.edge_curved:
            graph.es["edge_curved"] = options.edge_curved
        else:
            graph.es["edge_curved"] = 0.0
        if not options.file_layout:
            if options.use_boundary_edges:
                boundary_edges = get_boundary_edges(graph)
                for edge in boundary_edges:
                    graph.es[edge.index]['opacity'] = "rgba(0,0,0,0.1)"
            compute_layout(graph, options, boundary_edges=boundary_edges)

        visual_style = {}

        visual_style['edge_label'] = None
        visual_style['edge_color'] = graph.es['opacity']
        visual_style["edge_curved"] = graph.es["edge_curved"]
        visual_style['edge_width'] = graph.es['weight']

        visual_style['vertex_frame_width'] = graph.vs['vertex_frame_width']
        visual_style['vertex_shape'] = graph.vs['vertex_shape']
        visual_style['vertex_size'] = graph.vs['weight']
        visual_style['vertex_color'] = graph.vs['vertex_color']
        visual_style['vertex_frame_color'] = graph.vs['vertex_frame_color']

        visual_style['layout'] = graph.vs['layout']

        visual_style['bbox'] = options.bbox
        visual_style['margin'] = options.margin
        visual_style['edge_order_by'] = ('weight', 'asc')

        if options.save_pdf:
            igraph.plot(graph, options.output + '.pdf', **visual_style)
            if options.pdf_rotate:
                helper.rotate_pdf(options.output)
            if options.img_trim:
                command = 'pdfcrop ' + options.output + '.pdf ' + options.output + '.pdf'
                os.system(command)

        if options.save_png:
            igraph.plot(graph, options.output + '.png', **visual_style)
            if options.img_trim:
                command = 'convert ' + options.output + '.png -trim ' + options.output + '.png'
                os.system(command)

        if options.show_plot:
            igraph.plot(graph, **visual_style)
