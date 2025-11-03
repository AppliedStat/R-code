import numpy as np
import time
import random
from  functools import reduce
from scipy.spatial.distance import cdist
import operator
from scipy.spatial import distance
import matplotlib.pyplot as plt


#随机产生
#data = np.random.rand(2000, 2) * 100

#TSP-Lib数据集
file_path = r'C:\Users\USER\Desktop\GA解p-centera\pr439.txt' #，pr439 rat575 rat783 pr1002,rl1323
data = np.loadtxt(file_path, delimiter=" ")

def calculate_column_difference(matrix):
    # 将二维数组转换为NumPy数组
    arr = np.array(matrix)
    # 计算每一列的最大值和最小值
    max_values = np.amax(arr, axis=0)
    min_values = np.amin(arr, axis=0)
    # 计算差值
    differences = max_values - min_values
    return differences

# 计算每一列最大值和最小值的差值
column_diffs = calculate_column_difference(data)

class ProjectorStack:
    """
    Stack of points that are shifted / projected to put first one at origin.
    """
    def __init__(self, vec):
        self.vs = np.array(vec)
        
    def push(self, v):
        if len(self.vs) == 0:
            self.vs = np.array([v])
        else:
            self.vs = np.append(self.vs, [v], axis=0)
        return self
    
    def pop(self):
        if len(self.vs) > 0:
            ret, self.vs = self.vs[-1], self.vs[:-1]
            return ret
    
    def __mul__(self, v):
        s = np.zeros(len(v))
        for vi in self.vs:
            s = s + vi * np.dot(vi, v)
        return s
    
class GaertnerBoundary:
    """
        GärtnerBoundary

    See the passage regarding M_B in Section 4 of Gärtner's paper.
    """
    def __init__(self, pts):
        self.projector = ProjectorStack([])
        self.centers, self.square_radii = np.array([]), np.array([])
        self.empty_center = np.array([np.NaN for _ in pts[0]])


def push_if_stable(bound, pt):
    if len(bound.centers) == 0:
        bound.square_radii = np.append(bound.square_radii, 0.0)
        bound.centers = np.array([pt])
        return True
    q0, center = bound.centers[0], bound.centers[-1]
    C, r2  = center - q0, bound.square_radii[-1]
    Qm, M = pt - q0, bound.projector
    Qm_bar = M * Qm
    residue, e = Qm - Qm_bar, sqdist(Qm, C) - r2
    z, tol = 2 * sqnorm(residue), np.finfo(float).eps * max(r2, 1.0)
    isstable = np.abs(z) > tol
    if isstable:
        center_new  = center + (e / z) * residue
        r2new = r2 + (e * e) / (2 * z)
        bound.projector.push(residue / np.linalg.norm(residue))
        bound.centers = np.append(bound.centers, np.array([center_new]), axis=0)
        bound.square_radii = np.append(bound.square_radii, r2new)
    return isstable

def pop(bound):
    n = len(bound.centers)
    bound.centers = bound.centers[:-1]
    bound.square_radii = bound.square_radii[:-1]
    if n >= 2:
        bound.projector.pop()
    return bound


class NSphere:
    def __init__(self, c, sqr):
        self.center = np.array(c)
        self.sqradius = sqr

def isinside(pt, nsphere, atol=1e-6, rtol=0.0):
    r2, R2 = sqdist(pt, nsphere.center), nsphere.sqradius
    return r2 <= R2 or np.isclose(r2, R2, atol=atol**2,rtol=rtol**2)

def allinside(pts, nsphere, atol=1e-6, rtol=0.0):
    for p in pts:
        if not isinside(p, nsphere, atol, rtol):
            return False
    return True

def move_to_front(pts, i):
    pt = pts[i]
    for j in range(len(pts)):
        pts[j], pt = pt, np.array(pts[j])
        if j == i:
            break
    return pts

def dist(p1, p2):
    return np.linalg.norm(p1 - p2)

def sqdist(p1, p2):
    return sqnorm(p1 - p2)

def sqnorm(p):
    return np.sum(np.array([x * x for x in p]))

def ismaxlength(bound):
    len(bound.centers) == len(bound.empty_center) + 1

def makeNSphere(bound):
    if len(bound.centers) == 0: 
        return NSphere(bound.empty_center, 0.0)
    return NSphere(bound.centers[-1], bound.square_radii[-1])

def _welzl(pts, pos, bdry):
    support_count, nsphere = 0, makeNSphere(bdry)
    if ismaxlength(bdry):
        return nsphere, 0
    for i in range(pos):
        if not isinside(pts[i], nsphere):
            isstable = push_if_stable(bdry, pts[i])
            if isstable:
                nsphere, s = _welzl(pts, i, bdry)
                pop(bdry)
                move_to_front(pts, i)
                support_count = s + 1
    return nsphere, support_count

def find_max_excess(nsphere, pts, k1):
    err_max, k_max = -np.Inf, k1 - 1
    for (k, pt) in enumerate(pts[k_max:]):
        err = sqdist(pt, nsphere.center) - nsphere.sqradius
        if  err > err_max:
            err_max, k_max = err, k + k1
    return err_max, k_max - 1

#最小圆的圆心是平面任意一点
def welzl(points, maxiterations=2000):
    pts, eps = np.array(points, copy=True), np.finfo(float).eps
    bdry, t = GaertnerBoundary(pts), 1
    nsphere, s = _welzl(pts, t, bdry)
    for i in range(maxiterations):
        e, k = find_max_excess(nsphere, pts, t + 1)
        if e <= eps:
            break
        pt = pts[k]
        push_if_stable(bdry, pt)
        nsphere_new, s_new = _welzl(pts, s, bdry)
        pop(bdry)
        move_to_front(pts, k)
        nsphere = nsphere_new
        t, s = s + 1, s_new + 1
    return nsphere

def euclidean_distance(a, b):
    distance = np.linalg.norm(a-b, ord=2)
    return distance
"""
#最小圆的圆心在数据点上
def welzl(data):
    n = len(data)

    # If there are no points, return None
    if n == 0:
        return None

    # If there's only one point, the circle is centered on that point with radius 0
    if n == 1:
        return data[0], 0

    # Iterate over each point and calculate the circle centered on that point
    smallest_circle = None
    for i in range(n):
        center = data[i]
        radius = 0
        # Check if all other points are within the circle
        valid_circle = True
        for j in range(n):
            if j != i:
                if euclidean_distance(data[j], center) > radius:
                    # Update the radius if the current point is outside the circle
                    radius = euclidean_distance(data[j], center)

        # Update the smallest circle if necessary
        if valid_circle and (smallest_circle is None or radius < smallest_circle[1]):
            smallest_circle = (center, radius)

    return smallest_circle
"""

def sum_euclidean_distance(X, Y):
    return np.sqrt(np.sum(np.power(X - Y, 2), axis=1))

"""
#1.随机选择m个点
def init_populations(points, m, num_populations):
    populations = np.zeros((num_populations, m*2))
    
    for i in range(num_populations):
        center_indices = random.sample(range(len(points)), m)
        centers = points[center_indices].flatten()
        populations[i] = centers
    return populations
"""
def init_solutions(data, m, num_populations):
    #populations = np.zeros((num_populations, m))
    population_centers = []
    pop_radii = []
    fitnesses = np.zeros(num_populations)
    
    for j in range(num_populations):
        max_distances = np.zeros(m)
        # Initialize centroids list with first random centroid
        centers = [data[np.random.randint(data.shape[0])]]
        while len(centers) < m:
            dists = cdist(data, centers)
            min_dists = np.min(dists, axis=1)
            next_center = data[np.argmax(min_dists)]
            centers.append(next_center)
        # 计算每个数据点到中心点的距离
        distances = np.linalg.norm(data[:, np.newaxis, :] - centers, axis=2)
        # 找到每个数据点所属的簇索引（即距离最近的中心点索引）
        cluster_indices = np.argmin(distances, axis=1)
        # 计算每个簇中的最远数据点与对应中心的距离
        for i in range(len(data)):
            distance_to_center = distances[i, cluster_indices[i]]
            if distance_to_center > max_distances[cluster_indices[i]]:
                max_distances[cluster_indices[i]] = distance_to_center
        
        population_centers.append(np.array(centers))#.flatten()
        pop_radii.append(max_distances)
        fitnesses[j] = -max(max_distances)
    return np.array(population_centers), np.array(pop_radii), fitnesses

#2.#每次设置最远点(第一个点随机选取)
def init_populations(data, m, num_populations):
    #populations = np.zeros((num_populations, m))
    populations = []
    for i in range(num_populations):
        # Initialize centroids list with first random centroid
        centers = [data[np.random.randint(data.shape[0])]]
        while len(centers) < m:
            dists = cdist(data, centers)
            min_dists = np.min(dists, axis=1)
            next_center = data[np.argmax(min_dists)]
            centers.append(next_center)
        
        populations.append(np.array(centers))#.flatten()
    return np.array(populations)

def is_right_turn(p1, p2, p3):
    # Check if the three points make a right turn
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]) < 0

def melkman_algorithm(points):
    # Sort the points based on x-coordinate
    #sorted_points = sorted(points)
    sorted_points = points[np.lexsort(points[:,::-1].T)]
    # Initialize the output lists for upper and lower hulls
    upper_hull = []
    lower_hull = []

    # Build the upper hull
    for point in sorted_points:
        while len(upper_hull) >= 2 and is_right_turn(upper_hull[-2], upper_hull[-1], point):
            upper_hull.pop()
        upper_hull.append(point)

    # Build the lower hull
    for point in reversed(sorted_points):
        while len(lower_hull) >= 2 and is_right_turn(lower_hull[-2], lower_hull[-1], point):
            lower_hull.pop()
        lower_hull.append(point)

    # Combine the upper and lower hulls
    hull = upper_hull + lower_hull[1:-1]

    return np.array(hull)

def circumcircle(p_1, p_2, p_3):
    """
    :return:  x0 and y0 is center of a circle, r is radius of a circle
    """
    x1, y1 = p_1
    x2, y2 = p_2
    x3, y3 = p_3
    x0 = 1/2*((x1**2+y1**2)*(y3-y2)+(x2**2+y2**2)*(y1-y3)+(x3**2+y3**2)*(y2-y1))/(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))
    y0 = 1/2*((x1**2+y1**2)*(x3-x2)+(x2**2+y2**2)*(x1-x3)+(x3**2+y3**2)*(x2-x1))/(y1*(x3-x2)+y2*(x1-x3)+y3*(x2-x1))
    center = np.array([x0,y0])
    radius = distance.euclidean(center, p_1)
    return center, radius

def get_minimum_enclosing_circle(points):
    # Step 1
    distances = distance.cdist(points, points)
    max_distance_idx = np.unravel_index(np.argmax(distances), distances.shape)
    p1 = points[max_distance_idx[0]]
    p2 = points[max_distance_idx[1]]
    radius = distance.euclidean(p1, p2) / 2
    center = (p1 + p2) / 2

    # Check if all points are covered
    if np.all(np.round(distance.cdist([center], points),5) <= np.round(radius,5)):
        return center

    # Step 2 and Step 3
    else:
        center_distances = distance.cdist([center], points)
        max_cendis_idx = np.unravel_index(np.argmax(center_distances), center_distances.shape)
        p3 = points[max_cendis_idx[1]]
        center, radius = circumcircle(p1, p2, p3)

        # Check if all points are covered
        #print(distance.cdist([center], points),radius)
        if np.all(np.round(distance.cdist([center], points),5) <= np.round(radius,5)):
            return center
        
        else:
            center_distances = distance.cdist([center], points)
            max_cendis_idx = np.unravel_index(np.argmax(center_distances), center_distances.shape)
            p4 = points[max_cendis_idx[1]]
            
            #print(p1,p2,p3,p4)
            # Step 4
            d1 = distance.euclidean(p1, p3)
            d2 = distance.euclidean(p1, p4)
            d3 = distance.euclidean(p2, p3)
            d4 = distance.euclidean(p2, p4)
            
            # Choose the appropriate point to omit
            if min(d1, d2, d3, d4) == d1 or min(d1, d2, d3, d4) == d2:
                center, radius = circumcircle(p2, p3, p4)
                return center
            else:
                center, radius = circumcircle(p1, p3, p4)
                return center

def cluster_data_simple(cluster_center):
    global num_clusters, data
    cluster_index = np.argmin(np.linalg.norm(data[:, np.newaxis] - cluster_center, axis=2), axis=1)
    new_cluster_center = []
    for idx in range(num_clusters):
        clust_idx = np.where(cluster_index == idx)[0]
        if len(clust_idx) == 0:
            continue
        elif len(clust_idx) == 1:
            center = data[clust_idx][0]
        elif len(clust_idx) == 2:
            #center, radius = two_point_circle(data[clust_idx])
            center = (data[clust_idx][0]+data[clust_idx][1])/2
        else:
            one_cluster_data = data[clust_idx]
            #圆心在平面任意位置
            if len(one_cluster_data) > 100:
                one_cluster_data = melkman_algorithm(one_cluster_data)
                center = get_minimum_enclosing_circle(one_cluster_data)
            else:
                center = get_minimum_enclosing_circle(one_cluster_data)
        new_cluster_center.append(center)
    return np.array(new_cluster_center)

#给定中心坐标，给出每个点对应坐标的索引，最大距离等信息
def cluster_data_centers(cluster_center):
    global num_clusters, data 
    # 初始化每个簇的最远距离为0
    max_distances = np.zeros(len(cluster_center))
    # 计算每个数据点到中心点的距离
    distances = np.linalg.norm(data[:, np.newaxis, :] - cluster_center, axis=2)
    # 找到每个数据点所属的簇索引（即距离最近的中心点索引）
    cluster_indices = np.argmin(distances, axis=1)
    # 计算每个簇中的最远数据点与对应中心的距离
    for i in range(len(data)):
        distance_to_center = distances[i, cluster_indices[i]]
        if distance_to_center > max_distances[cluster_indices[i]]:
            max_distances[cluster_indices[i]] = distance_to_center

    return cluster_indices, max_distances, -max(max_distances)

def cluster_data(solution):
    global num_clusters, feature_vector_length, data
    """
    cluster_centers = []
    for clust_idx in range(num_clusters):
        cluster_centers.append(solution[feature_vector_length*clust_idx:feature_vector_length*(clust_idx+1)])
    cluster_centers = np.array(cluster_centers)
    """
    index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - solution, axis=2), axis=1)
    new_centers = []
    radii = []
    for idx in range(num_clusters):
        clust_idx = np.where(index_sequence == idx)[0]
        if len(clust_idx) == 0:
            # for the first reset_generation
            #center = 0
            #radius = 0
            # for the second reset_generation
            #radius = 0
            #radii.append(radius)
            # for the thrid reset_generation
            continue
        elif len(clust_idx) == 1:
            center = data[clust_idx][0]
            radius = 0
        elif len(clust_idx) == 2:
            #center, radius = two_point_circle(data[clust_idx])
            center = (data[clust_idx][0]+data[clust_idx][1])/2
            radius = np.linalg.norm(center - data[clust_idx][0])
        else:
            one_cluster_data = data[clust_idx]
            #圆心在平面任意位置
            nsphere = welzl(one_cluster_data)
            center = nsphere.center
            radius = np.sqrt(nsphere.sqradius)
        #圆心在一个点上
        #center, radius = nsphere
        ##clusters.append(np.where(solution == clust_idx)[0])
        new_centers.append(center)
        radii.append(radius)
    fitness = -max(radii)
    
    return np.array(new_centers), np.array(radii), fitness

def reset_func(solutions):
    global num_parents_mating
    #num_solutions = len(solutions)
    sol_centers = []
    max_dists = []
    fitnesses = np.zeros(num_parents_mating)
    for i in range(num_parents_mating):
        center, radii, fitness = cluster_data(solutions[i])
        sol_centers.append(center)
        max_dists.append(radii)
        fitnesses[i] = fitness
    #fitness = 1.0 / (maxmin_dist + 0.00000001)
    return np.array(sol_centers), np.array(max_dists), fitnesses
    

#稳态选择，对应第二个__init__，效果更好
def parent_selection_func(populations, fitness, num_parents_mating):
    global num_genes
    fitness_sorted = sorted(range(len(fitness)), key=lambda k: fitness[k])
    fitness_sorted.reverse()

    part_parents = np.empty((num_parents_mating, num_genes))

    for parent_num in range(num_parents_mating):
        """"
        _, cluster_centers = cluster_data(ga_instance.population[fitness_sorted[parent_num], :].copy(), fitness_sorted[parent_num])
        
        assignments = []
        for point in data:
            distances = np.linalg.norm(point - cluster_centers, axis=1)
            closest_center_index = np.argmin(distances)
            assignments.append(closest_center_index)
        parents[parent_num, :] = np.array(assignments)
        """
        part_parents[parent_num, :] = populations[fitness_sorted[parent_num], :].copy()

    return part_parents, fitness_sorted[:num_parents_mating]

#自己设计的精英插入选择
def parent_selection_func_1(populations, fitness, num_parents):
    global num_genes
    fitness_sorted = sorted(range(len(fitness)), key=lambda k: fitness[k])
    fitness_sorted.reverse()
    #print(fitness)

    num_parents_mating = int(num_parents/3)
    rows_to_remove = np.random.choice(populations.shape[0], num_parents_mating, replace=False)
    parents_1 = np.delete(populations, rows_to_remove, axis=0)
    fitness_1 = np.delete(fitness, rows_to_remove)
    parents_2 = np.empty((num_parents_mating, populations.shape[1]))
    for parent_num in range(num_parents_mating):
        parents_2[parent_num, :] = populations[fitness_sorted[parent_num], :].copy()
    parents = np.concatenate((parents_1, parents_2), axis=0)

    fitness = np.append(fitness_1, fitness[fitness_sorted[:num_parents_mating]])
    #print(fitness)
    fitness_sorted = sorted(range(len(fitness)), key=lambda k: fitness[k])
    fitness_sorted.reverse()
    
    #求fitness但从第二轮开始出错，不知道为什么
    """
    fitness = []
    for index, individuality in enumerate(parents):
        one_fitness = fitness_func(parents, individuality, index)
        fitness.append(one_fitness)

    fitness_sorted = sorted(range(len(fitness)), key=lambda k: fitness[k])
    fitness_sorted.reverse()
    """
    return parents, fitness_sorted

#每个基因被变异的概率相同
def mutation_func_1(offspring):
    global num_genes, data
    Pm = 0.09
    for chromosome_idx in range(offspring.shape[0]):
        for gene_idx in range(num_genes):
            prob = np.random.uniform(0.0, 1.0)
            if prob < Pm:
                # 生成一个在[-0.1, 0.1]范围内的随机小数
                mutation_amount = np.random.uniform(-0.1, 0.1)
                if gene_idx % 2 == 0:
                    mutation_amount *= column_diffs[0]
                elif gene_idx % 2 == 1:
                    mutation_amount *= column_diffs[1]
                offspring[chromosome_idx, gene_idx] += mutation_amount
    return offspring

#半径越大，基因被变异的概率越大
def mutation_func(offsprings, circles_dist):
    global num_parents_mating, num_clusters, data
    #Pm = 2
    for chromosome_idx in range(num_parents_mating):
        #prob = np.random.uniform(0.0, 1.0)
        #if prob < Pm:
        #_, cluster_centers, radii = cluster_data(offsprings[chromosome_idx])
        #gap = 1/circles_dist[chromosome_idx]
        gap = circles_dist[chromosome_idx]
        
        """
        gene_idx = random.randint(0, num_clusters-1)
        #1.随机簇的索引，距离越大，选中的概率越高
        Pr = gap / sum(gap)
        gene_idx = np.random.choice(range(num_clusters), p=Pr.ravel())
        """
        #1.最大簇的索引
        gene_idx = np.argmax(gap)
        #gene_idx = np.argmin(gap)
        #gene_idx = random.randint(0, num_clusters-1)
        dimensions = data.shape[1]
        
        mutation_number = random.randint(0, dimensions-1)
            # 生成一个在[-0.1, 0.1]范围内的随机小数
        mutation_amount = np.random.uniform(-0.1, 0.1)
        mutation_amount *= column_diffs[mutation_number]
        offsprings[chromosome_idx][gene_idx][mutation_number] += mutation_amount
        """
        for i in range(dimensions):
            # 生成一个在[-0.1, 0.1]范围内的随机小数
            mutation_amount = np.random.uniform(-0.1, 0.1)
            mutation_amount *= column_diffs[i]
            offsprings[chromosome_idx][gene_idx][i] += mutation_amount
        """
    return offsprings

def reset_generation_1(populations):
    global num_parents_mating, num_clusters, data
    circles_radii = []
    
    reset_populs = []
    reset_fitnesses = []
    for chromo_idx in range(populations.shape[0]):
        while True:
            sol_center, max_distances, fitness = cluster_data(populations[chromo_idx])
            index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - populations[chromo_idx], axis=2), axis=1)
                    
            unique_elements, counts = np.unique(index_sequence, return_counts=True)
            
            if len(unique_elements) == num_clusters:
                reset_populs.append(sol_center)
                circles_radii.append(max_distances)
                reset_fitnesses.append(fitness)
                break
            else:
                num_largest = num_clusters - len(unique_elements)
                largest_indices = np.argsort(max_distances)[-num_largest:][::-1]
                
                # 查找哪些中心点没有分配到任何需求点
                unassigned_centers = set(range(num_clusters))
                assigned_centers = set(unique_elements)
                unassigned_centers -= assigned_centers
                unassigned_centers = list(unassigned_centers)
                m = len(unassigned_centers)
                for i in range(m):
                    clust_id = np.where(index_sequence == largest_indices[i])[0]
                    coordinates = data[clust_id]
                    random_center = coordinates[np.random.choice(coordinates.shape[0])]
                    #incom_centers = np.vstack((incom_centers, random_center))
                    populations[chromo_idx][unassigned_centers[i]] = random_center
                
                #populations[chromo_idx],  max_distances, reset_fitness[chromo_idx] = cluster_data(populations[chromo_idx])
                #circles_radii.append(max_distances)
    return np.array(reset_populs), np.array(circles_radii), np.array(reset_fitnesses)

def reset_generation_2(populations):
    global num_parents_mating, num_clusters, data
    circles_radii = []
    
    reset_populs = []
    reset_fitnesses = []
    for chromo_idx in range(populations.shape[0]):
        sol_center = populations[chromo_idx]
        while True:
            index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - sol_center, axis=2), axis=1)
            new_sol_center, max_distances, fitness = cluster_data(sol_center)
            if len(new_sol_center) == num_clusters:
                reset_populs.append(new_sol_center)
                circles_radii.append(max_distances)
                reset_fitnesses.append(fitness)
                break
            else:
                num_largest = num_clusters - len(new_sol_center)
                largest_indices = np.argsort(max_distances)[-num_largest:][::-1]
                
                #index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - new_sol_center, axis=2), axis=1)      
                for i in range(num_largest):
                    clust_id = np.where(index_sequence == largest_indices[i])[0]
                    coordinates = data[clust_id]
                    random_center = coordinates[np.random.choice(coordinates.shape[0])]
                    sol_center = np.vstack((new_sol_center, random_center))
                    #sol_center = np.append(sol_center, [random_center], axis=0)
                #populations[chromo_idx] = sol_center

    return np.array(reset_populs), np.array(circles_radii), np.array(reset_fitnesses)

def reset_generation(populations):
    global num_parents_mating, num_clusters, data
    circles_radii = []
    
    reset_populs = []
    reset_fitnesses = []
    for chromo_idx in range(populations.shape[0]):
        sol_center = populations[chromo_idx]
        while True:
            new_sol_center = cluster_data_simple(sol_center)
            index_sequence, max_distances, fitness = cluster_data_centers(new_sol_center)
            if len(new_sol_center) == num_clusters:
                reset_populs.append(new_sol_center)
                circles_radii.append(max_distances)
                reset_fitnesses.append(fitness)
                break
            else:
                num_largest = num_clusters - len(new_sol_center)
                #index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - new_sol_center, axis=2), axis=1)      
                for i in range(num_largest):
                    dists = cdist(data, new_sol_center)
                    min_dists = np.min(dists, axis=1)
                    next_center = data[np.argmax(min_dists)]
                    sol_center = np.vstack((new_sol_center, next_center))
                    #sol_center = np.append(sol_center, [random_center], axis=0)
                #populations[chromo_idx] = sol_center

    return np.array(reset_populs), np.array(circles_radii), np.array(reset_fitnesses)

def reset_generation_备份(populations):
    global num_parents_mating, num_clusters, data
    circles_radii = []
    
    reset_populs = []
    reset_fitnesses = []
    for chromo_idx in range(populations.shape[0]):
        sol_center = populations[chromo_idx]
        while True:
            new_sol_center = cluster_data_simple(sol_center)
            index_sequence, max_distances, fitness = cluster_data_centers(new_sol_center)
            if len(new_sol_center) == num_clusters:
                reset_populs.append(new_sol_center)
                circles_radii.append(max_distances)
                reset_fitnesses.append(fitness)
                break
            else:
                num_largest = num_clusters - len(new_sol_center)
                largest_indices = np.argsort(max_distances)[-num_largest:][::-1]
                
                #index_sequence = np.argmin(np.linalg.norm(data[:, np.newaxis] - new_sol_center, axis=2), axis=1)      
                for i in range(num_largest):
                    clust_id = np.where(index_sequence == largest_indices[i])[0]
                    coordinates = data[clust_id]
                    random_center = coordinates[np.random.choice(coordinates.shape[0])]
                    sol_center = np.vstack((new_sol_center, random_center))
                    #sol_center = np.append(sol_center, [random_center], axis=0)
                #populations[chromo_idx] = sol_center

    return np.array(reset_populs), np.array(circles_radii), np.array(reset_fitnesses)


#最后统计各种结果
def statics(value, global_min):

    min_value = min(value)
    max_value = max(value)
    value = np.array(value)
    dev = (value-global_min)/global_min
    min_dev = min(dev)
    max_dev = max(dev)
    dev1 = np.sum(dev<0.05)
    dev2 = np.sum(dev<0.1)
    return min_value, max_value, min_dev, max_dev, dev1, dev2

"""
#将全部重新reset_generation
if __name__ == "__main__":
    
    num_clusters = 5
    feature_vector_length = data.shape[1]
    num_genes = num_clusters * feature_vector_length
    num_iterations=50
    num_populations=10 #种群数目
    
    start_time = time.time()
    populations = init_populations(data, num_clusters, num_populations) #构造初始种群
    set_fitness = fitness_func(populations)
    for i in range(num_iterations):   
        select_parents, _ = parent_selection_func(populations, set_fitness, num_populations)
        #new_offspring = mutation_func(select_parents)
        populations = mutation_func(select_parents)
        #稳态将最后几个父代替换掉
        #remain_num = num_populations - num_parents_mating
        #if remain_num > 0:
            #reser_parents, reser_fitness_index = parent_selection_func(populations, set_fitness, remain_num)
        #populations = np.concatenate((reser_parents, new_offspring), axis=0)

        populations = reset_generation(populations)
        set_fitness = fitness_func(populations)
        best_solution = populations[np.argmax(set_fitness)]
        best_fitness = max(set_fitness)

"""
if __name__ == "__main__":
    
    num_clusters = 10
    feature_vector_length = data.shape[1]
    num_genes = num_clusters * feature_vector_length
    #num_iterations=50
    num_populations=31 #种群数目
    num_parents_mating=10
    
    best_results = []
    each_time = []
    
    for i in range(100):
        
        start_time = time.time()
        #1.开始构造最小圆
        #populations = init_populations(data, num_clusters, num_populations)
        #parents, parent_radii, set_fitness = reset_func(populations)
        #1.开始构造最小圆，而是直接最远随机点
        parents, parent_radii, set_fitness = init_solutions(data, num_clusters, num_populations)
        #print(i, max(set_fitness))
        #init_fitness = float("inf")
        #count = 0
        #for i in range(num_iterations):
        #while True:
        for j in range(100):
            #select_parents, _ = parent_selection_func(populations, set_fitness, num_parents_mating)
     
            mating_indices = np.argsort(-set_fitness)[:num_parents_mating]
            select_parents = parents[mating_indices]
            select_radii = parent_radii[mating_indices]
            
            new_offsprings = mutation_func(select_parents, select_radii)
            reset_parents, reset_radii, reset_fitness = reset_generation(new_offsprings)
            #稳态将最后几个父代替换掉
            remain_num = num_populations - num_parents_mating
            if remain_num > 0:
                #reser_parents, reser_fitness_index = parent_selection_func(populations, set_fitness, remain_num)
                #reser_fitness = set_fitness[np.isin(np.arange(len(set_fitness)), reser_fitness_index)]
            
                reser = np.argsort(-set_fitness)[:remain_num]
                reser_parents = parents[reser]
                reser_radii = parent_radii[reser]
                reser_fitness = set_fitness[reser]
                #print(len(reset_indices))
                #new_offspring = reset_generation(new_offspring)
                #select_fitness = fitness_func(new_offspring)
                parents = np.concatenate((reser_parents, reset_parents), axis=0)
                parent_radii = np.concatenate((reser_radii, reset_radii), axis=0)
                set_fitness = np.concatenate([reser_fitness, reset_fitness])
            else:
                parents = reset_parents
                parent_radii = reset_radii
                set_fitness = reset_fitness
          
        end_time = time.time()
        run_time = end_time - start_time
        each_time.append(run_time)
        best_fitness = -max(set_fitness)
        print(i, best_fitness)
        best_results.append(best_fitness)
        best_solution = parents[np.argmax(set_fitness)]
    
    best, worst, min_dev, max_dev, dev1, dev2 = statics(best_results, 706.145)
    print(best, worst, min_dev, max_dev, dev1, dev2)
    print(each_time)
    print(np.mean(best_results), np.mean(each_time))

    min_dists, cluster_centers, _ = cluster_data(best_solution)
    
    #print(max(min_dists), cluster_centers)
    
    cluster_centers = []
    all_clusters_dists = []
    clusters = []
    
    for clust_idx in range(num_clusters):
        
        cluster_centers.append(best_solution[feature_vector_length*clust_idx:feature_vector_length*(clust_idx+1)])
        cluster_center_dists = sum_euclidean_distance(data, cluster_centers[clust_idx])
        all_clusters_dists.append(np.array(cluster_center_dists))
    
    cluster_centers = np.array(cluster_centers)
    all_clusters_dists = np.array(all_clusters_dists)
    cluster_indices = np.argmin(all_clusters_dists, axis=0)
    
    for clust_idx in range(num_clusters):
        clusters.append(np.where(cluster_indices == clust_idx)[0])
    
    ax = plt.subplot(111)
    for cluster_idx in range(num_clusters):
        cluster_x = data[clusters[cluster_idx], 0]
        cluster_y = data[clusters[cluster_idx], 1]
        plt.scatter(cluster_x, cluster_y)
        plt.scatter(cluster_centers[cluster_idx, 0], cluster_centers[cluster_idx, 1], linewidths=5)
        cir1 = plt.Circle((cluster_centers[cluster_idx, 0], cluster_centers[cluster_idx, 1]), max(min_dists), color="black",fill=False)
        ax.add_patch(cir1)
    plt.title("Clustering using PyGAD")
    plt.xlim(-5.5,20.5)
    plt.ylim(-5.5,20.5)
    plt.xticks([-5,0,5,10,15,20])
    plt.yticks([-5,0,5,10,15,20])
    ax.set_aspect('equal', adjustable='box')
    plt.show()
