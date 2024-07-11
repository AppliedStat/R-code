#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:14:33 2020

@author: OS
"""

# from __future__ import print_function
# Linux: 
# pip install weightedstats 
# pip install matplotlib 
# pip install sklearn
# pip install scikit-learn
# pip install numpy 

import numpy as np
import random
import pandas as pd
import weightedstats as ws
from noisify.recipes import machine_error,human_error
from functools import reduce
import matplotlib.pyplot as plt
import scipy.optimize as optim
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA

def creat_weight(n):
    #random.seed(2)
    nums = [random.uniform(0,1) for x in range(0,n)]
    summary = reduce(lambda x,y: x+y, nums)
    norm = [x/summary for x in nums]
    weight_set = np.array(norm)
    return weight_set

#Lp distance between two points
def get_distance(a, b, p):
    distance = np.linalg.norm(a-b, ord=p)
    return distance

def sum_distance(theta, x_y):
    #for i in range(len(x_y)):
        #dist = np.sum(w[i] * np.linalg.norm(x_y[i]-theta[i], ord=k))
    #x_y = np.column_stack((x,y))
    num = x_y.shape[0]
    dist = 0.0
    for i in range(num):
        distmat = np.linalg.norm(x_y[i]-theta[i], ord=2)
        #distmat = w[i]*distmat
        dist = dist + distmat
    return dist

#objective function
def weighed_distance(theta, x_y, w):
    #for i in range(len(x_y)):
        #dist = np.sum(w[i] * np.linalg.norm(x_y[i]-theta[i], ord=k))
    #x_y = np.column_stack((x,y))
    num = x_y.shape[0]
    dist = 0.0
    for i in range(num):
        distmat = np.linalg.norm(x_y[i]-theta, ord=2)
        distmat = w[i]*distmat
        dist = dist + distmat
    return dist

#两个数组元素相乘求和
def sumAndMul2List(list1, list2):

    result = np.sum([a*b for a,b in zip(list1,list2)])

    return result

#Optimizing the center points with BFGS
def k_weighted_means(samples, clusters, k, cutoff):
    #samples_xy = samples[:,0:2].tolist()
    #clusters = random.sample(samples_xy, k)
    #clusters = np.array(clusters)
    #clusters = np.array([[-3,0],[0,-0.5],[3,0.5]])
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            smallest_distance = get_distance(sample[0:2], clusters[0],2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1],2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])  
        biggest_shift = 0.0
        for j in range(k):
            #print(len(lists_xy[j]))
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            clusters_x = ws.weighted_mean(lists_xy[j][:,0], lists_w[j])
            clusters_y = ws.weighted_mean(lists_xy[j][:,1], lists_w[j])

            spare_clusters = np.append(spare_clusters, [[clusters_x, clusters_y]], axis=0)
            shift = get_distance(clusters[j], [clusters_x, clusters_y],2)
            biggest_shift = max(biggest_shift, shift)
        #print(biggest_shift)
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
            #print(lists_w)
    
    return clusters, lists_xy, lists_w

#Sorting to determine the weighted median
def w_median(points, w_value):
    Z = zip(points, w_value)
    Z = sorted(Z)
    points_sorted, w_value_sorted = zip(*Z)
    points_sorted = np.array(points_sorted)
    w_value_sorted = np.array(w_value_sorted)
    return points_sorted, w_value_sorted

#Sum of the previous subarrays
def part_sum(len_k, w_sorted):
    p_sum = 0
    for i in range(len_k):
        p_sum += w_sorted[i]
        if p_sum > 0.50000000000000000 or p_sum == 0.50000000000000000:
            return i
        else:
            continue
     
#Optimizing the center points with Manhattan
def k_weighted_medians(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            #print(len(lists_xy[j]))
            #part_sum = 0
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            clusters_x = ws.weighted_median(lists_xy[j][:,0], lists_w[j])
            clusters_y = ws.weighted_median(lists_xy[j][:,1], lists_w[j])
                              
            spare_clusters = np.append(spare_clusters, [[clusters_x, clusters_y]], axis=0)
            shift = get_distance(clusters[j], [clusters_x, clusters_y], 2)
            biggest_shift = max(biggest_shift, shift)
        #print(biggest_shift)
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
            #print(lists_w)
    return clusters, lists_xy, lists_w

#Optimizing the center points with BFGS
def bdfs_l2(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            #print(sample)
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            j_clusters = clusters[j]
            j_clusters.resize((1,2))
            #lists_xy[j][:,0].resize((:,1))
            len_lists = len(lists_xy[j])
            w_lists = lists_w[j]
            w_lists.resize((len_lists,1))
            #print(j_clusters)
            #print(lists_xy[j])
            #print(w_lists)
            #result = optim.fmin_l_bfgs_b(weighed_distance, clusters[j], args=(lists_xy[j][:,0], lists_xy[j][:,1], lists_w[j]))
            #init_theta = np.random.normal(size=(1, 2))
            #print(init_theta)
            #dis = weighed_distance(init_theta, lists_xy[j], w_lists)
            #print(lists_xy[j])
            result = optim.fmin_l_bfgs_b(weighed_distance, j_clusters, args=(lists_xy[j], w_lists), approx_grad=True)
            
            spare_clusters = np.append(spare_clusters, [[result[0][0], result[0][1]]], axis=0)
            shift = get_distance(clusters[j], result[0], 2)
            biggest_shift = max(biggest_shift, shift)
            
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
    return clusters, lists_xy, lists_w

#K-weighted-L2-median,用weiszfeld算法实现
def weight_weiszfeld_algorithm(points, normalized_weights, epsilon=1e-6, max_iterations=100):
    
    # Normalize the weights
    #normalized_weights = weights / np.sum(weights)
    # Initialize the center
    center = np.average(points, axis=0, weights=normalized_weights)
    
    for _ in range(max_iterations):
        # Calculate the distances from each point to the current center
        distances = np.linalg.norm(points - center, axis=1)
        
        # Handle zero distances
        distances = np.where(distances == 0, epsilon, distances)
        
        # Update the center based on weighted average
        new_center = np.sum(points * (normalized_weights / distances)[:, np.newaxis], axis=0) / np.sum(normalized_weights / distances)
        
        # Check for convergence
        if np.linalg.norm(new_center - center) < epsilon:
            break
        
        center = new_center
    
    return center

def k_weights_weiszfelds(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            #print(sample)
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            center = weight_weiszfeld_algorithm(lists_xy[j],lists_w[j])
            
            spare_clusters = np.append(spare_clusters, [center], axis=0)
            shift = get_distance(clusters[j], center,2)
            
            biggest_shift = max(biggest_shift, shift)
            
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
    return clusters, lists_xy, lists_w

def HL_estimator(criter, array_point):
    n = len(array_point)
    k_hl = []
    if criter == 'hl_1':
        for i in range(n-1):
            for j in range(i+1, n):
                med_value = (array_point[i]+array_point[j])/2
                k_hl.append(med_value)
    elif criter == 'hl_2':
        for i in range(n):
            for j in range(i, n):
                med_value = (array_point[i]+array_point[j])/2
                k_hl.append(med_value)
    elif criter == 'hl_3':
        for i in range(n):
            for j in range(n):
                med_value = (array_point[i]+array_point[j])/2
                k_hl.append(med_value)
                
    k_hl = np.array(k_hl)
    return k_hl

#Optimizing the center points with Manhattan
def k_weighted_hl1s(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            #print(len(lists_xy[j]))
            #part_sum = 0
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            lists_xhl = HL_estimator('hl_1', lists_xy[j][:,0])
            lists_yhl = HL_estimator('hl_1', lists_xy[j][:,1])
            lists_whl = HL_estimator('hl_1', lists_w[j])
            clusters_x = ws.weighted_median(lists_xhl, lists_whl)
            clusters_y = ws.weighted_median(lists_yhl, lists_whl)

            spare_clusters = np.append(spare_clusters, [[clusters_x, clusters_y]], axis=0)
            shift = get_distance(clusters[j], [clusters_x, clusters_y], 2)
            biggest_shift = max(biggest_shift, shift)
        #print(biggest_shift)
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
            #print(lists_w)
    return clusters, lists_xy, lists_w

def k_weighted_hl2s(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            #print(len(lists_xy[j]))
            #part_sum = 0
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            lists_xhl = HL_estimator('hl_2', lists_xy[j][:,0])
            lists_yhl = HL_estimator('hl_2', lists_xy[j][:,1])
            lists_whl = HL_estimator('hl_2', lists_w[j])
            clusters_x = ws.weighted_median(lists_xhl, lists_whl)
            clusters_y = ws.weighted_median(lists_yhl, lists_whl)
            spare_clusters = np.append(spare_clusters, [[clusters_x, clusters_y]], axis=0)
            shift = get_distance(clusters[j], [clusters_x, clusters_y], 2)
            biggest_shift = max(biggest_shift, shift)
        #print(biggest_shift)
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
            #print(lists_w)
    return clusters, lists_xy, lists_w

def k_weighted_hl3s(samples, clusters, k, cutoff):
    n_loop = 0
    while True:
        #print(clusters)
        spare_clusters = np.empty(shape=[0, 2])
        lists_xy = [[] for _ in range(k)]
        lists_w = [[] for _ in range(k)]
        for sample in samples:
            smallest_distance = get_distance(sample[0:2], clusters[0], 2)
            cluster_index = 0
            
            for i in range(k-1):
                distance = get_distance(sample[0:2], clusters[i+1], 2)
                if distance < smallest_distance:
                    smallest_distance = distance
                    cluster_index = i + 1
                    
            lists_xy[cluster_index].append(sample[0:2])
            lists_w[cluster_index].append(sample[2])
            
        biggest_shift = 0.0
        for j in range(k):
            #print(len(lists_xy[j]))
            #part_sum = 0
            lists_xy[j] = np.array(lists_xy[j])
            lists_w[j] = np.array(lists_w[j])
            lists_xhl = HL_estimator('hl_3', lists_xy[j][:,0])
            lists_yhl = HL_estimator('hl_3', lists_xy[j][:,1])
            lists_whl = HL_estimator('hl_3', lists_w[j])
            clusters_x = ws.weighted_median(lists_xhl, lists_whl)
            clusters_y = ws.weighted_median(lists_yhl, lists_whl)
            spare_clusters = np.append(spare_clusters, [[clusters_x, clusters_y]], axis=0)
            shift = get_distance(clusters[j], [clusters_x, clusters_y], 2)
            biggest_shift = max(biggest_shift, shift)
        #print(biggest_shift)
        if biggest_shift < cutoff:
            #print("第{}次迭代后，聚类稳定。".format(n_loop))
            break
        else:
            clusters = spare_clusters
            #print(clusters)
            n_loop += 1
            #print(lists_w)
    return clusters, lists_xy, lists_w

def rela_eff(x,a1,b1,a2,b2):
    sum_x11 = 0
    sum_x12 = 0
    sum_x11_x12 = 0
    sum_x21 = 0
    sum_x22 = 0
    sum_x21_x22 = 0
    for i in range(len(x)):
        x11 = x[i][0][0]-a1
        x12 = x[i][0][1]-b1
        sum_x11 += pow(x11,2)
        sum_x12 += pow(x12,2)
        sum_x11_x12 += x11*x12
        
        x21 = x[i][1][0]-a2
        x22 = x[i][1][1]-b2
        sum_x21 += pow(x21,2)
        sum_x22 += pow(x22,2)
        sum_x21_x22 += x21*x22
    return (sum_x11*sum_x12-pow(sum_x11_x12,2))/pow(len(x),2),(sum_x21*sum_x22-pow(sum_x21_x22,2))/pow(len(x),2)

def re_assign(samples, clusters, metric):
    k = len(clusters)
    lists_xy = [[] for _ in range(k)]
    for sample in samples:
        smallest_distance = get_distance(sample[0:2], clusters[0], metric)
        cluster_index = 0
        for i in range(k-1):
            distance = get_distance(sample[0:2], clusters[i+1], metric)
            if distance < smallest_distance:
                smallest_distance = distance
                cluster_index = i + 1
        lists_xy[cluster_index].append(sample[0:2])
    return lists_xy

def multisource_weber_objective(points, centers, weights):

    total_weighted_distance = 0
    for point, weight in zip(points, weights):
        distances = np.linalg.norm(point - centers, axis=1)
        min_distance = np.min(distances)
        total_weighted_distance += weight * min_distance
    return total_weighted_distance


if __name__ == '__main__':

    all_clusters_means = []
    all_clusters_medians = []
    all_clusters_bfgs = []
    all_clusters_hl1 = []
    all_clusters_hl2 = []
    all_clusters_hl3 = []
    
    all_clusters_means_c1 = []
    all_clusters_means_c2 = []
    all_clusters_medians_c1 = []
    all_clusters_medians_c2 = []
    all_clusters_bfgs_c1 = []
    all_clusters_bfgs_c2 = []
    all_clusters_hl1_c1 = []
    all_clusters_hl1_c2 = []
    all_clusters_hl2_c1 = []
    all_clusters_hl2_c2 = []
    all_clusters_hl3_c1 = []
    all_clusters_hl3_c2 = []
    
    mean_distance = []
    median_distance = []
    l2_median_distance = []
    hl_distance = []
    # 随机生成三组二元正态分布随机数
    #np.random.seed(1314)
    for i in range(1):
        #np.random.seed(2)

        filename = "SJC1.txt"     # <-------------------------------------------
        point_array = np.loadtxt(filename)
        
        #Laplace
        x1, y1 = np.random.laplace(1, 1, 100).T, np.random.laplace(1, 1, 100).T
        x2, y2 = np.random.laplace(3, 1, 100).T, np.random.laplace(3, 1, 100).T

        #x_out = np.array([0,3,x2[0],10,30,50,70])
        #y_out = np.array([0,3,y2[0],10,30,50,70])
        #for j,k in zip(x_out, y_out):

        # 绘制三组数据的散点图
        x = np.hstack((x1, x2))
        y = np.hstack((y1, y2))
        #print(x,y)

        iris = load_iris()
        data = iris.data    # <-------------------------------------
        #data[:] = data[:] + np.random.random()-0.5
        
        #scaler = StandardScaler()
        #data = scaler.fit_transform(wine)
        
        data_copy = data.copy()
        
        #对数据加入污染
        #data[0:10]= data[0:10]*2
        #data[:] = machine_noise(data[:])
        
        pca = PCA(n_components=2)#实例化
        pca = pca.fit(data)#拟合模型
        data = pca.transform(data)#获取新矩阵

        pca = PCA(n_components=2)#实例化
        pca = pca.fit(data_copy)#拟合模型
        data_copy = pca.transform(data_copy)#获取新矩阵
        #w_set = creat_weight(len(data))
        w_set = [1 for i in range(len(data))]

        x_y_w = np.column_stack((data,w_set))
        x_y_w_copy = np.column_stack((data_copy,w_set))

        k=3

        clusters = np.array([[-3,0],[0,-0.5],[3,0.5]])
        #clusters = np.array([[-400,0],[0,0],[600,0]])
        #clusters = np.array([[-4,1],[1,-1],[5,0]])

        no_conta = []
        conta = []
        
        n_clusters_means, coor_xy_means, coor_w_means= k_weighted_means(x_y_w, clusters, k, 0.00000001)
        no_conta.append(n_clusters_means)
        #print(n_clusters_means)
        n_clusters_means_copy, _, _= k_weighted_means(x_y_w_copy, clusters, k, 0.00000001)
        conta.append(n_clusters_means_copy)
        #print(n_clusters_means_copy)
        dist_mean = sum_distance(n_clusters_means, n_clusters_means_copy)
        mean_distance.append(dist_mean)
       
        n_clusters_medians, coor_xy_medians, coor_w_medians= k_weighted_medians(x_y_w, clusters, k, 0.00000001)
        no_conta.append(n_clusters_medians)
        #print(n_clusters_medians)
        n_clusters_medians_copy, _, _= k_weighted_medians(x_y_w_copy, clusters, k, 0.00000001)
        conta.append(n_clusters_medians_copy)
        
        dist_median = sum_distance(n_clusters_medians, n_clusters_medians_copy)
        median_distance.append(dist_median)
 
        n_clusters_bfgs, coor_xy_bfgs, coor_w_bfgs= k_weights_weiszfelds(x_y_w, clusters, k, 0.00000001)
        no_conta.append(n_clusters_bfgs)
        #print(n_clusters_bfgs)
        n_clusters_bfgs_copy, _, _= k_weights_weiszfelds(x_y_w_copy, clusters, k, 0.00000001)
        conta.append(n_clusters_bfgs_copy)
        
        dist_l2_median = sum_distance(n_clusters_bfgs, n_clusters_bfgs_copy)
        l2_median_distance.append(dist_l2_median)

        n_clusters_hl1, coor_xy_hl1, coor_w_hl1= k_weighted_hl1s(x_y_w, clusters, k, 0.00000001)
        no_conta.append(n_clusters_hl1)
        #print(n_clusters_hl1)
        n_clusters_hl1_copy, _, _ = k_weighted_hl1s(x_y_w_copy, clusters, k, 0.00000001)
        conta.append(n_clusters_hl1_copy)
        
        dist_hl = sum_distance(n_clusters_hl1, n_clusters_hl1_copy)
        hl_distance.append(dist_hl)
                   
        #n_clusters_hl2, coor_xy_hl2, coor_w_hl2= hod_mann_hl2(x_y_w, clusters, k, 0.00000001)
        #print(n_clusters_hl2)
        #n_clusters_hl2_copy, _, _= hod_mann_hl2(x_y_w_copy, clusters, k, 0.00000001)

        #n_clusters_hl3, coor_xy_hl3, coor_w_hl3= hod_mann_hl3(x_y_w, clusters, k, 0.00000001)
        #print(n_clusters_hl3)
        #n_clusters_hl3_copy, _, _= hod_mann_hl3(x_y_w_copy, clusters, k, 0.00000001)

        #print(coor_xy)
        all_clusters_means.append(n_clusters_means)
        all_clusters_medians.append(n_clusters_medians)
        all_clusters_bfgs.append(n_clusters_bfgs)
        all_clusters_hl1.append(n_clusters_hl1)

    print(np.mean(mean_distance))
    print(np.mean(median_distance))
    print(np.mean(l2_median_distance))
    print(np.mean(hl_distance))
    
    coor_means = re_assign(data_copy, n_clusters_means,2)
    print(len(coor_means[0]),len(coor_means[1]),len(coor_means[2]))
    coor_medians = re_assign(data_copy, n_clusters_medians,2)
    print(len(coor_medians[0]),len(coor_medians[1]),len(coor_medians[2]))
    coor_bfgs = re_assign(data_copy, n_clusters_bfgs,2)
    print(len(coor_bfgs[0]),len(coor_bfgs[1]),len(coor_bfgs[2]))
    #coor_bfgs_ = re_assign(data_copy, n_clusters_bfgs_,2)
    #print(len(coor_bfgs_[0]),len(coor_bfgs_[1]),len(coor_bfgs_[2]))
    coor_hl1 = re_assign(data_copy, n_clusters_hl1,2)
    print(len(coor_hl1[0]),len(coor_hl1[1]),len(coor_hl1[2]))
    #coor_hl2 = re_assign(data_copy, n_clusters_hl2,2)
    #coor_hl3 = re_assign(data_copy, n_clusters_hl3,2)
    
    #画出分类图
    plt.style.use('ggplot')
    color = ["red","blue","orange"]
    for i in [0,1]:
        coor_means[i] = np.array(coor_means[i])
        plt.scatter(coor_means[i][:,0],
                    coor_means[i][:,1],
                    s=10,
                    c=color[i],
                    alpha = 0.5)#透明度
    #plt.xlim(-3.5,4)
    plt.show()
    #plt.scatter(coor_xy_means[0][:,0],coor_xy_means[0][:,1],c="red")
    #plt.scatter(coor_xy_means[1][:,0],coor_xy_means[1][:,1],c = "black")
    #plt.scatter(coor_xy_means[2][:,0],coor_xy_means[2][:,1],c="orange")
    #plt.legend()#显示图例
    #plt.title("PCA of IRIS dataset")#显示标题
    
    for i in [0,1,2]:
        coor_medians[i] = np.array(coor_medians[i])
        plt.scatter(coor_medians[i][:,0],
                    coor_medians[i][:,1],
                    c=color[i],
                    s=10,
                    alpha = 0.5)#透明度
    #plt.xlim(-3.5,4)
    plt.show()
    
    for i in [0,1,2]:
        coor_bfgs[i] = np.array(coor_bfgs[i])
        plt.scatter(coor_bfgs[i][:,0],
                    coor_bfgs[i][:,1],
                    c=color[i],
                    s=10,
                    alpha = 0.5)#透明度
    #plt.xlim(-3.5,4)
    plt.show()
    
    for i in [0,1,2]:
        coor_hl1[i] = np.array(coor_hl1[i])
        plt.scatter(coor_hl1[i][:,0],
                    coor_hl1[i][:,1],
                    c=color[i],
                    s=10,
                    alpha = 0.5)#透明度
    #plt.xlim(-3.5,4)
    plt.show()    ## Figure 
    
#plt.grid(axis="y")

