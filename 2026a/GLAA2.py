# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 00:15:00 2023

@author: GAO
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 21:59:48 2023

@author: GAO
"""

import numpy as np
import time
import random
from scipy.spatial import distance
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt


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

#最小圆的圆心在平面的任意位置
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
"""
def euclidean_distance(a, b):
    distance = np.linalg.norm(a-b, ord=2)
    return distance

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
def set_euclidean_distance(X, Y):
    return np.sqrt(np.sum(np.power(X - Y, 2), axis=1))

"""
#1.随机选择m个不重复点
def init_solutions(data, m, num_populations):
    #populations = np.zeros((num_populations, m))
    population_centers = []
    pop_cluster_indices = []
    pop_radii = []
    fitnesses = np.zeros(num_populations)
    
    for j in range(num_populations):
        max_distances = np.zeros(m)
        # Initialize centroids list with first random centroid
        center_indices = random.sample(range(len(data)), m)
        centers = data[center_indices]
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
        pop_cluster_indices.append(cluster_indices)
        pop_radii.append(max_distances)
        fitnesses[j] = -max(max_distances)
        
    return np.array(pop_cluster_indices), np.array(population_centers), np.array(pop_radii), fitnesses

"""
#2.#每次设置最远点(第一个点随机选取)
def init_solutions(data, m, num_populations):
    #populations = np.zeros((num_populations, m))
    population_centers = []
    pop_cluster_indices = []
    pop_radii = []
    fitnesses = np.zeros(num_populations)
    first_points = random.sample(range(data.shape[0]), num_populations)
    
    for j in range(num_populations):
        max_distances = np.zeros(m)
        # Initialize centroids list with first random centroid
        #centers = [data[np.random.randint(data.shape[0])]]
        centers = [data[first_points[j]]]
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
        pop_cluster_indices.append(cluster_indices)
        pop_radii.append(max_distances)
        fitnesses[j] = -max(max_distances)
    return np.array(pop_cluster_indices), np.array(population_centers), np.array(pop_radii), fitnesses

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

#只通过数据点的索引求中心坐标
def cluster_data_simple(cluster_index):
    global num_clusters, data
    #cluster_index = np.argmin(np.linalg.norm(data[:, np.newaxis] - circle_center, axis=2), axis=1)
    cluster_center = []
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
            if len(one_cluster_data) > 100:
                one_cluster_data = melkman_algorithm(one_cluster_data)
                center = get_minimum_enclosing_circle(one_cluster_data)
            else:
                center = get_minimum_enclosing_circle(one_cluster_data)
            #圆心在平面任意位置
            #nsphere = welzl(one_cluster_data)
            #center = nsphere.center
        cluster_center.append(center)
    return np.array(cluster_center)

#通过数据点的索引求中心坐标，最大距离以及下一个迭代的数据点索引等信息，连接reset_func
def cluster_data(cluster_index):
    global num_clusters, data
    cluster_center = []
    radii = []
    for idx in range(num_clusters):
        clust_idx = np.where(cluster_index == idx)[0]
        if len(clust_idx) == 0:
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
            one_cluster_data =  melkman_algorithm(one_cluster_data)
            #圆心在平面任意位置
            nsphere = welzl(one_cluster_data)
            center = nsphere.center
            radius = np.sqrt(nsphere.sqradius)
        cluster_center.append(center)
        radii.append(radius)
    cluster_center = np.array(cluster_center)
    #new_cluster_index = np.argmin(np.linalg.norm(data[:, np.newaxis] - cluster_center, axis=2), axis=1)
    fitness = -max(radii)
    return cluster_center, radii, fitness

def reset_func(cluster_indexes):
    global num_parents_mating
    #new_cluster_indexes = []
    cluster_centers = []
    max_dists = []
    fitnesses = np.zeros(num_parents_mating)
    for i in range(num_parents_mating):
        cluster_center, radii, fitness = cluster_data(cluster_indexes[i])
        #new_cluster_indexes.append(new_cluster_index)
        cluster_centers.append(cluster_center)
        max_dists.append(radii)
        fitnesses[i] = fitness
    return np.array(cluster_centers), np.array(max_dists), fitnesses

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

#稳态选择,对应第二个__init__,效果更好。我们此次代码在主函数中以更简略的方式直接实现
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


#修改变异过程，根据Genetic K-Means Algorithm 这篇论文
def mutation_func_1(offspring, centers):
    global num_genes, data
    Pm = 0.05
    #cm = 2
    for chromosome_idx in range(offspring.shape[0]):
        for gene_idx in range(num_genes):
            prob = np.random.uniform(0.0, 1.0)
            if prob < Pm:
                #centers, _, _ = cluster_data(offspring[chromosome_idx])
                gene_clusters_dists = set_euclidean_distance(data[gene_idx],centers)
                if np.count_nonzero(gene_clusters_dists) == len(gene_clusters_dists):
                    #max_dist = max(gene_clusters_dists)
                    #gap = cm*max_dist - gene_clusters_dists
                    #Pr = gap / sum(gap)
                    gap = 1/gene_clusters_dists
                    Pr = gap / sum(gap)
                    new_gene = np.random.choice(range(num_clusters), p=Pr.ravel())
                    offspring[chromosome_idx, gene_idx] = new_gene       
        #random_gene_idx = np.random.choice(range(offspring.shape[1]))
        #offspring[chromosome_idx, random_gene_idx] += np.random.random()
    return offspring

def mutation_func(offsprings, centers, circles_dist):
    global num_clusters, data
    for chromosome_idx in range(offsprings.shape[0]):
        #选择发生变异的簇
        which_cluster_idx = np.argmax(circles_dist[chromosome_idx])
        #which_cluster_idx = np.argmin(circles_dist[chromosome_idx])
        #which_cluster_idx = random.randint(0, num_clusters-1)
        #选择簇中的某个点
        clust_idx = np.where(offsprings[chromosome_idx] == which_cluster_idx)[0]
        #if len(clust_idx) != 1:   #只在测试最小簇和随机簇时使用
        clus_data = data[clust_idx]
        inverse_distances = set_euclidean_distance(clus_data, centers[chromosome_idx][which_cluster_idx])
        #print(clust_idx, clus_data, centers[chromosome_idx][which_cluster_idx], inverse_distances)
        #inverse_distances = 1.0 / (1.0 + clusdata_center_dists)  # 距离越小，概率越大
        #inverse_distances = clusdata_center_dists
        Pd = inverse_distances / np.sum(inverse_distances)
        
        new_idx = np.random.choice(range(len(inverse_distances)), p=Pd.ravel())
        
        chosen_idx = clust_idx[new_idx]
        
        #离哪个索引越近，变异成哪个索引的概率越大
        gene_clusters_dists = set_euclidean_distance(data[chosen_idx], centers[chromosome_idx])
        new_gene_idx = np.argsort(gene_clusters_dists)[1]
        offsprings[chromosome_idx][chosen_idx] = new_gene_idx

    return offsprings

#修改变异过程，根据Genetic K-Means Algorithm 这篇论文，但是选择的某点的索引不能变异为原来的索引
def mutation_func_2(offspring):
    global num_genes, data
    Pm = 0.0001
    cm = 2
    for chromosome_idx in range(offspring.shape[0]):
        for gene_idx in range(num_genes):
            prob = np.random.uniform(0.0, 1.0)
            if prob < Pm:
                _, centers, _ = cluster_data(offspring[chromosome_idx])
                centers = centers[~np.any(centers == centers[offspring[chromosome_idx, gene_idx]], axis=1)]
                gene_clusters_dists = set_euclidean_distance(data[gene_idx],centers)
                if np.count_nonzero(gene_clusters_dists) == len(gene_clusters_dists):
                    #print(gene_clusters_dists)
                    max_dist = max(gene_clusters_dists)
                    gap = cm*max_dist - gene_clusters_dists
                    #gap = gene_clusters_dists
                    Pr = gap / sum(gap)
                    #print(np.delete(range(num_clusters), offspring[chromosome_idx, gene_idx]),Pr)
                    new_gene = np.random.choice(np.delete(range(num_clusters), offspring[chromosome_idx, gene_idx]), p=Pr.ravel())
                    offspring[chromosome_idx, gene_idx] = new_gene       
        #random_gene_idx = np.random.choice(range(offspring.shape[1]))
        #offspring[chromosome_idx, random_gene_idx] += np.random.random()
    return offspring


#修改变异过程，选取每个簇内的某个点，距离中心越远，被选中概率越大
def mutation_func_3(offspring):
    global num_clusters, data, num_genes
    Pm = 0.02
    cm = 1
    for chromosome_idx in range(offspring.shape[0]):
        for num_clusters_idx in range(num_clusters):
            prob = np.random.uniform(0.0, 1.0)
            if prob < Pm:
                #选择发生变异的点
                clust_idx = np.where(offspring[chromosome_idx] == num_clusters_idx)[0]
                if len(clust_idx) == 0 or len(clust_idx) == 1:
                    continue
                clus_data = data[clust_idx]
                nsphere = welzl(clus_data)
                center, radius = nsphere
                data_center_dists = set_euclidean_distance(clus_data,center)
                Pr = data_center_dists / sum(data_center_dists)
                new_idx = np.random.choice(range(len(data_center_dists)), p=Pr.ravel())
                chosen_idx = clust_idx[new_idx]
                
                #变异成哪个索引
                _, centers, _ = cluster_data(offspring[chromosome_idx])
                #删除之前的圆心
                centers = centers[~np.any(centers == center, axis=1)]
                gene_clusters_dists = set_euclidean_distance(data[chosen_idx],centers)
                if np.count_nonzero(gene_clusters_dists) == len(gene_clusters_dists):
                    max_dist = max(gene_clusters_dists)
                    gap = cm*max_dist - gene_clusters_dists
                    #gap = gene_clusters_dists
                    Pr = gap / sum(gap)
                    new_gene_idx = np.random.choice(np.delete(range(num_clusters), num_clusters_idx), p=Pr.ravel())
                    offspring[chromosome_idx, chosen_idx] = new_gene_idx 
    return offspring

#修改变异过程，选取单点，争取选大半径圆内的点且距离中心越远，该点被选中概率越大
def mutation_func_4(offsprings, centers, circles_dist):
    global num_clusters, data
    for chromosome_idx in range(offsprings.shape[0]):
        #选择发生变异的簇
        Pc = circles_dist[chromosome_idx] / sum(circles_dist[chromosome_idx])
        which_cluster_idx = np.random.choice(range(num_clusters), p=Pc.ravel())
        #选择簇中的某个点
        clust_idx = np.where(offsprings[chromosome_idx] == which_cluster_idx)[0]
        if len(clust_idx) == 0 or len(clust_idx) == 1:
            continue
        clus_data = data[clust_idx]
        data_center_dists = set_euclidean_distance(clus_data,centers[chromosome_idx][which_cluster_idx])
        Pd = data_center_dists / sum(data_center_dists)
        new_idx = np.random.choice(range(len(data_center_dists)), p=Pd.ravel())
        chosen_idx = clust_idx[new_idx]
        
        #离哪个索引越近，变异成哪个索引的概率越大
        #删除之前的圆心
        #centers = centers[~np.any(centers == centers[which_cluster_idx], axis=1)]
        gene_clusters_dists = set_euclidean_distance(data[chosen_idx], centers[chromosome_idx])
        if np.count_nonzero(gene_clusters_dists) == len(gene_clusters_dists):
            gap = 1/gene_clusters_dists
            Pr = gap / sum(gap)
            #new_gene_idx = np.random.choice(np.delete(range(num_clusters), which_cluster_idx), p=Pr.ravel())
            new_gene_idx = np.random.choice(range(num_clusters), p=Pr.ravel())
            offsprings[chromosome_idx][chosen_idx] = new_gene_idx

    return offsprings

def reset_generation(populations):
    global num_clusters, data
    reset_populs = []
    circles_radii = []
    reset_fitnesses = []
    reset_centers = []
    for chromo_idx in range(populations.shape[0]):
        
        #cluster_indexes = populations[chromo_idx]
        incom_centers = cluster_data_simple(populations[chromo_idx])
        #new_cluster_indexes, incom_centers, radiuses, fitness = cluster_data(cluster_indexes)
        while True: 
            #new_cluster_indexes = np.argmin(np.linalg.norm(data[:, np.newaxis] - incom_centers, axis=2), axis=1)
            new_cluster_indexes, radiuses, fitness = cluster_data_centers(incom_centers)
            unique_elements, _ = np.unique(new_cluster_indexes, return_counts=True)
            if len(unique_elements) == num_clusters:
                reset_populs.append(new_cluster_indexes)
                reset_centers.append(incom_centers)
                circles_radii.append(radiuses)
                reset_fitnesses.append(fitness)
                break
            else:
                #reset_popul.append(populations[chromo_idx])
                num_largest = num_clusters - len(unique_elements)
                # 查找哪些中心点没有分配到任何需求点
                unassigned_centers = list(set(range(num_clusters)) - set(new_cluster_indexes))

                for i in range(num_largest):
                    dists = cdist(data, incom_centers)
                    min_dists = np.min(dists, axis=1)
                    next_center = data[np.argmax(min_dists)]
                    if len(incom_centers) == num_clusters:
                        incom_centers[unassigned_centers[i]] = next_center
                    else:
                        incom_centers = np.vstack((incom_centers, next_center))

    return np.array(reset_populs), np.array(reset_centers), np.array(circles_radii), np.array(reset_fitnesses)


def reset_generation_备份(populations):
    global num_clusters, data
    reset_populs = []
    circles_radii = []
    reset_fitnesses = []
    reset_centers = []
    for chromo_idx in range(populations.shape[0]):
        
        cluster_indexes = populations[chromo_idx]
        incom_centers = cluster_data_simple(cluster_indexes)
        #new_cluster_indexes, incom_centers, radiuses, fitness = cluster_data(cluster_indexes)
        while True: 
            #new_cluster_indexes = np.argmin(np.linalg.norm(data[:, np.newaxis] - incom_centers, axis=2), axis=1)
            new_cluster_indexes, radiuses, fitness = cluster_data_centers(incom_centers)
            unique_elements, _ = np.unique(new_cluster_indexes, return_counts=True)
            if len(unique_elements) == num_clusters:
                reset_populs.append(new_cluster_indexes)
                reset_centers.append(incom_centers)
                circles_radii.append(radiuses)
                reset_fitnesses.append(fitness)
                break
            else:
                #reset_popul.append(populations[chromo_idx])
                num_largest = num_clusters - len(unique_elements)
                largest_indices = np.argsort(radiuses)[-num_largest:][::-1]
                
                # 查找哪些中心点没有分配到任何需求点
                unassigned_centers = list(set(range(num_clusters)) - set(new_cluster_indexes))

                for i in range(num_largest):
                    clust_id = np.where(new_cluster_indexes == largest_indices[i])[0]
                    coordinates = data[clust_id]
                    random_center = coordinates[np.random.choice(coordinates.shape[0])]
                    if len(incom_centers) == num_clusters:
                        incom_centers[unassigned_centers[i]] = random_center
                    else:
                        incom_centers = np.vstack((incom_centers, random_center))

    return np.array(reset_populs), np.array(reset_centers), np.array(circles_radii), np.array(reset_fitnesses)

def locate_allocate(populations):
    global num_clusters, data
    reset_populs = []
    circles_radii = []
    reset_fitnesses = []
    reset_centers = []
    for chromo_idx in range(populations.shape[0]):
        iterate_radii = []
        cluster_indexes = populations[chromo_idx]
        while True:
            incom_centers = cluster_data_simple(cluster_indexes)
            new_cluster_indexes, radiuses, fitness = cluster_data_centers(incom_centers)
            if any(np.array_equal(radiuses, arr) for arr in iterate_radii):
                break
            else:
                iterate_radii.append(radiuses)
        
        reset_populs.append(new_cluster_indexes)
        reset_centers.append(incom_centers)
        circles_radii.append(radiuses)
        reset_fitnesses.append(fitness)

    return np.array(reset_populs), np.array(reset_centers), np.array(circles_radii), np.array(reset_fitnesses)

def reset_generation_1(populations):
    global num_clusters, data
    reset_populs = []
    circles_radii = []
    reset_fitnesses = []
    reset_centers = []
    for chromo_idx in range(populations.shape[0]):
        
        cluster_indexes = populations[chromo_idx]
        incom_centers = cluster_data_simple(cluster_indexes)
        while True:
            new_cluster_indexes = np.argmin(np.linalg.norm(data[:, np.newaxis] - incom_centers, axis=2), axis=1)
            update_centers, radiuses, fitness = cluster_data(new_cluster_indexes)
            if len(update_centers) == num_clusters:
                reset_populs.append(new_cluster_indexes)
                reset_centers.append(update_centers)
                circles_radii.append(radiuses)
                reset_fitnesses.append(fitness)
                break
            else:
                #reset_popul.append(populations[chromo_idx])
                num_largest = num_clusters - len(incom_centers)
                largest_indices = np.argsort(radiuses)[-num_largest:][::-1]

                for i in range(num_largest):
                    clust_id = np.where(new_cluster_indexes == largest_indices[i])[0]
                    coordinates = data[clust_id]
                    random_center = coordinates[np.random.choice(coordinates.shape[0])]
                    incom_centers = np.vstack((update_centers, random_center))
                    #incom_centers[unassigned_centers[i]] = random_center

    return np.array(reset_populs), np.array(reset_centers), np.array(circles_radii), np.array(reset_fitnesses)

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
    num_iterations=50
    num_populations=20 #种群数目
    num_parents_mating=5
    num_genes = len(data)
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

if __name__ == "__main__":
    
    num_clusters = 5
    num_iterations=5
    num_populations=3 #种群数目
    num_parents_mating=2
    num_genes = len(data)
    start_time = time.time()
    populations = init_populations(data, num_clusters, num_populations) #构造初始种群
    set_fitness = fitness_func(populations)
    select_parents, _ = parent_selection_func(populations, set_fitness, num_parents_mating)
    for i in range(num_iterations):
        new_offspring = mutation_func(select_parents)
        new_offspring = reset_generation(new_offspring)
        select_fitness = fitness_func(new_offspring)
        merge_populations = np.concatenate((populations, new_offspring), axis=0)
        merge_fitness = np.concatenate([set_fitness, select_fitness])
        print(merge_fitness)
        max_indices = np.argsort(-merge_fitness)[:num_populations]
        populations = merge_populations[max_indices]
        set_fitness = merge_fitness[max_indices]
        print(set_fitness)
        best_solution = populations[np.argmax(set_fitness)]
        best_solution_1 = populations[:1]
        #print(best_solution == best_solution_1)
        select_parents = populations[:num_parents_mating]
        
        
#只将变异的重新reset_generation
if __name__ == "__main__":
    
    num_clusters = 5
    num_iterations=5
    num_populations=3 #种群数目
    num_parents_mating=2
    num_genes = len(data)
    start_time = time.time()
    populations = init_populations(data, num_clusters, num_populations) #构造初始种群
    set_fitness = fitness_func(populations)
    for i in range(num_iterations):   
        select_parents, _ = parent_selection_func(populations, set_fitness, num_parents_mating)
        new_offspring = mutation_func(select_parents)
        #稳态将最后几个父代替换掉
        remain_num = num_populations - num_parents_mating
        if remain_num > 0:
            reser_parents, reser_fitness_index = parent_selection_func(populations, set_fitness, remain_num)
        reser_fitness = set_fitness[np.isin(np.arange(len(set_fitness)), reser_fitness_index)]
        new_offspring = reset_generation(new_offspring)
        select_fitness = fitness_func(new_offspring)
        populations = np.concatenate((reser_parents, new_offspring), axis=0)
        set_fitness = np.concatenate([reser_fitness, select_fitness])
        #populations = reset_generation(populations)
        #set_fitness = fitness_func(populations)
        best_solution = populations[np.argmax(set_fitness)]
        #best_fitness = max(set_fitness)      
""" 
#只将变异的重新reset_generation
if __name__ == "__main__":
    
    #随机产生
    #data = np.random.rand(2000, 2) * 100
    
    #TSP-Lib数据集
    file_path = r'C:\Users\USER\Desktop\GA解p-centera\pr439.txt' #，pr439 rat575 rat783 pr1002  rl1323
    data = np.loadtxt(file_path, delimiter=" ")
    
    num_clusters = 100
    #num_iterations=40
    num_populations=31 #种群数目
    num_parents_mating=20
    num_genes = len(data)
    
    best_results = []
    each_time = []
    
    for i in range(100):
        #count = 1
        #init_fitness = float("inf")
        
        start_time = time.time()
        
        #parents_for_index = init_populations(data, num_clusters, num_populations)
        #parents_centers, parents_radii, set_fitnesses = reset_func(parents_for_index)
        parents_for_index, parents_centers, parents_radii, set_fitnesses = init_solutions(data, num_clusters, num_populations) #构造初始种群

        #parents_for_index, parents_centers, parents_radii, set_fitnesses = locate_allocate(parents_for_index)
        #while True:
        for j in range(100):
            #select_parents, _ = parent_selection_func(populations, set_fitness, num_parents_mating)
            
            mating_indices = np.argsort(-set_fitnesses)[:num_parents_mating]
            select_parents = parents_for_index[mating_indices]
            select_centers = parents_centers[mating_indices]
            select_radii = parents_radii[mating_indices]
            
            new_offsprings = mutation_func(select_parents, select_centers, select_radii)
            reset_parents, reset_centers, reset_radii, reset_fitness = reset_generation(new_offsprings)
            #reset_parents, reset_centers, reset_radii, reset_fitness = locate_allocate(reset_parents)
            
            #稳态将最后几个父代替换掉
            remain_num = num_populations - num_parents_mating
            if remain_num > 0:
                #reser_parents, reser_fitness_index = parent_selection_func(populations, set_fitness, remain_num)
                #reser_fitness = set_fitness[np.isin(np.arange(len(set_fitness)), reser_fitness_index)]
                reser = np.argsort(-set_fitnesses)[:remain_num]
                
                reser_parents = parents_for_index[reser]
                reser_centers = parents_centers[reser]
                reser_radii = parents_radii[reser]
                reser_fitness = set_fitnesses[reser]

                parents_for_index = np.concatenate((reser_parents, reset_parents), axis=0)
                parents_radii = np.concatenate((reser_radii, reset_radii), axis=0)
                parents_centers = np.concatenate((reser_centers, reset_centers), axis=0)
                set_fitnesses = np.concatenate([reser_fitness, reset_fitness])
                
            else:
                parents_for_index = reset_parents
                parents_centers = reset_centers
                parents_radii = reset_radii
                set_fitnesses = reset_fitness
            #populations = reset_generation(populations)
            #set_fitness = fitness_func(populations)
            #best_solution = populations[np.argmax(set_fitness)]
            #print(-max(set_fitnesses))
            best_fitness = -max(set_fitnesses)
            """
            if best_fitness == init_fitness:
                count += 1
                #if count == 2:
                    #num_parents_mating = int(num_populations/2)
                #elif count == 5:
                    #num_parents_mating = num_populations - 1
                if count == 2000:
                    best_solution = parents_for_index[np.argmax(set_fitnesses)]
                    break
            else:
                count = 1
                #num_parents_mating=2
                init_fitness = best_fitness
            """
        
        end_time = time.time()
        run_time = end_time - start_time
        #print("run_time:", run_time)
        each_time.append(run_time)
        best_fitness = -max(set_fitnesses)
        print(i, best_fitness) 
        best_results.append(best_fitness)
        best_solution = parents_for_index[np.argmax(set_fitnesses)]
        
    best, worst, min_dev, max_dev, dev1, dev2 = statics(best_results, 604.152)
    print(best, worst, min_dev, max_dev, dev1, dev2)
    print(each_time)
    print(np.mean(best_results), np.mean(each_time))
    
    cluster_centers, maxmin_dist, _ = cluster_data(best_solution)
    
    #print(maxmin_dist, cluster_centers)

    all_clusters_dists = []
    clusters = []
    
    for clust_idx in range(num_clusters):
        
        #cluster_centers.append(best_solution[feature_vector_length*clust_idx:feature_vector_length*(clust_idx+1)])
        cluster_center_dists = set_euclidean_distance(data, cluster_centers[clust_idx])
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
        cir1 = plt.Circle((cluster_centers[cluster_idx, 0], cluster_centers[cluster_idx, 1]), maxmin_dist, color="black",fill=False)
        ax.add_patch(cir1)
    plt.title("Clustering using PyGAD")
    plt.xlim(-5.5,20.5)
    plt.ylim(-5.5,20.5)
    plt.xticks([-5,0,5,10,15,20])
    plt.yticks([-5,0,5,10,15,20])
    ax.set_aspect('equal', adjustable='box')
    plt.show()
