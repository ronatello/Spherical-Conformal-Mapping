''' read the file line by line to understanding how the half edge load 
an off file and to save the half edge mesh as an obj file.

author:
    zhangsihao yang
'''
# print ('hello')

# data_path = 'tests/data/brain_python.off'
data_path = './tuette.off'

import mesh
from mesh import normalize, norm, dot, cross_product

# HalfedgeMesh
mesh = mesh.HalfedgeMesh(data_path)

# Returns a list of Vertex type (in order of file)--similarly for halfedges,
# and facets
mesh.vertices

# The number of facets in the mesh
# print(len(mesh.facets))

# Get the 10th halfedge
mesh.halfedges[10]

# Get the halfedge that starts at vertex 25 and ends at vertex 50
# mesh.get_halfedge(0, 1)

# print(mesh.vertices)
# for vertex in mesh.vertices:
#     print(vertex.get_vertex())

# print(mesh.facets)
# print('--------------------')
# for face in mesh.facets:
#     print(face.a, face.b, face.c)


def gauss_map(mesh):

    for vertex in mesh.vertices:

        he = vertex.halfedge
        face = he.facet
        normals = [face.get_normal()]
        he = he.opposite.prev

        while (he.facet != face):
            normals.append(he.facet.get_normal())
            he = he.opposite.prev

        vertex_normal = [sum(val) for val in zip(*normals)]
        vertex_normal = [val / len(normals) for val in vertex_normal]
        vertex_normal = normalize(vertex_normal)

        vertex.x = vertex_normal[0]
        vertex.y = vertex_normal[1]
        vertex.z = vertex_normal[2]

def tuette_energy(mesh):

    total_energy = 0

    for he in mesh.halfedges:
        v1 = he.vertex.get_vertex()
        v2 = he.opposite.vertex.get_vertex()

        vector_diff = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]

        energy = norm(vector_diff)
        energy = energy ** 2
        total_energy += energy

    return total_energy / 2

def harmonic_energy(mesh):

    total_energy = 0

    for he in mesh.halfedges:
        v1 = he.vertex.get_vertex()
        v2 = he.opposite.vertex.get_vertex()

        v3 = he.next.vertex.get_vertex()
        v4 = he.opposite.next.vertex.get_vertex()

        Ta = 0.5 * dot([v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]], [v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]]) * (1/(norm(cross_product([v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]], [v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]]))))
        Tb = 0.5 * dot([v1[0] - v4[0], v1[1] - v4[1], v1[2] - v4[2]], [v2[0] - v4[0], v2[1] - v4[1], v2[2] - v4[2]]) * (1/(norm(cross_product([v1[0] - v4[0], v1[1] - v4[1], v1[2] - v4[2]], [v2[0] - v4[0], v2[1] - v4[1], v2[2] - v4[2]]))))
        kuv = Ta + Tb

        vector_diff = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]

        energy = norm(vector_diff)
        energy = energy ** 2
        energy = kuv * energy

        total_energy += energy

    return total_energy / 2

def mobius_transformation(mesh):

    sphere_center = [0, 0, 0]

    for vertex in mesh.vertices:

        he = vertex.halfedge
        face = he.facet
        
        v1 = he.vertex.get_vertex()
        v2 = he.next.vertex.get_vertex()
        v3 = he.prev.vertex.get_vertex()

        v1_2 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]
        v1_3 = [v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]]

        area_of_triangle = 0.5 * norm(cross_product(v1_2, v1_3))
        area_element = 1/3 * area_of_triangle

        he = he.opposite.prev

        while (he.facet != face):

            v1 = he.vertex.get_vertex()
            v2 = he.next.vertex.get_vertex()
            v3 = he.prev.vertex.get_vertex()

            v1_2 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]
            v1_3 = [v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]]

            area_of_triangle = 0.5 * norm(cross_product(v1_2, v1_3))
            area_element += 1/3 * area_of_triangle

            he = he.opposite.prev

        sphere_center = [sphere_center[0] + (area_element * vertex.get_vertex()[0]), sphere_center[1] + (area_element * vertex.get_vertex()[1]), sphere_center[2] + (area_element * vertex.get_vertex()[2])]

    sphere_center = [val / len(mesh.vertices) for val in sphere_center]

    for vertex in mesh.vertices:

        vertex.x = vertex.x - sphere_center[0]
        vertex.y = vertex.y - sphere_center[1]
        vertex.z = vertex.z - sphere_center[2]

        vertex_norm = normalize(vertex.get_vertex())

        vertex.x = vertex_norm[0]
        vertex.y = vertex_norm[1]
        vertex.z = vertex_norm[2]

def absolute_derivative(mesh, step_size):

    update = []

    for vertex in mesh.vertices:

        he = vertex.halfedge
        head = he

        next_b_edge = he
        adjacent_vertices = []

        while next_b_edge.opposite:
            adjacent_vertices.append(next_b_edge.opposite.vertex.get_vertex())
            next_b_edge = next_b_edge.opposite.prev
            if (next_b_edge == head):
                break

        vertex_differences = [[vertex.get_vertex()[0] - val[0], vertex.get_vertex()[1] - val[1], vertex.get_vertex()[2] - val[2]] for val in adjacent_vertices]
        vertex_differential = [sum(val) for val in zip(*vertex_differences)]

        product_differential_normal = dot(vertex_differential, normalize(vertex.get_vertex()))
        norm_component = [product_differential_normal * val for val in normalize(vertex.get_vertex())]
        
        abs_der = [vertex_differential[0] - norm_component[0], vertex_differential[1] - norm_component[1], vertex_differential[2] - norm_component[2]]

        update.append([-val * step_size for val in abs_der])

    return update
            
def update_step(mesh, update):

    for i in range(len(mesh.vertices)):

        mesh.vertices[i].x = mesh.vertices[i].x + update[i][0]
        mesh.vertices[i].y = mesh.vertices[i].y + update[i][1]
        mesh.vertices[i].z = mesh.vertices[i].z + update[i][2]

        normed_vertex = normalize(mesh.vertices[i].get_vertex())

        mesh.vertices[i].x = normed_vertex[0]
        mesh.vertices[i].y = normed_vertex[1]
        mesh.vertices[i].z = normed_vertex[2]

def spherical_tuette_mapping(mesh, step_size, stop_value):

    initial_energy = tuette_energy(mesh)
    new_energy = 0

    while(initial_energy - new_energy > stop_value):
        update = absolute_derivative(mesh, step_size)
        update_step(mesh, update)
        if (new_energy != 0):
            initial_energy = new_energy
        new_energy = tuette_energy(mesh)

        print(new_energy)

def spherical_conformal_mapping(mesh, step_size, stop_value):

    initial_energy = tuette_energy(mesh)
    new_energy = 0

    print(initial_energy)

    while(initial_energy - new_energy > stop_value):
        update = absolute_derivative(mesh, step_size)
        update_step(mesh, update)
        mobius_transformation(mesh)
        if (new_energy != 0):
            initial_energy = new_energy
        new_energy = harmonic_energy(mesh)

        print(new_energy)

# to save the halfedge mesh you will need to following function.
def save_halfmesh_as_obj(mesh, file_name):
    with open(file_name, 'w') as open_file:
        for vertex in mesh.vertices:
            lv = vertex.get_vertex()
            open_file.write("v {} {} {} \n".format(lv[0], lv[1], lv[2]))

        for face in mesh.facets:
            open_file.write("f {} {} {}\n".format(face.a+1, face.b+1, face.c+1))

def custom_save_halfmesh_as_off(mesh, file_name):
    with open(file_name, 'w') as open_file:
        open_file.write('OFF\n')
        open_file.write(str(len(mesh.vertices)) + " " + str(len(mesh.facets)) + " 0" + "\n")
        for vertex in mesh.vertices:
            lv = vertex.get_vertex()
            open_file.write("{} {} {} \n".format(lv[0], lv[1], lv[2]))

        for face in mesh.facets:
            open_file.write("3 {} {} {}\n".format(face.a, face.b, face.c))

gauss_map(mesh)

save_halfmesh_as_obj(mesh, 'gauss.obj')
custom_save_halfmesh_as_off(mesh, 'gauss.off')

spherical_tuette_mapping(mesh, 0.01, 0.0001)

save_halfmesh_as_obj(mesh, 'tuette.obj')
custom_save_halfmesh_as_off(mesh, 'tuette.off')

spherical_conformal_mapping(mesh, 0.000001, 0.000001)

save_halfmesh_as_obj(mesh, 'conformal.obj')
custom_save_halfmesh_as_off(mesh, 'conformal.off')