import numpy as np
import timeit
from numba import cuda, jit
from tqdm import tqdm


USE_64 = True

if USE_64:
    bits = 64
    np_type = np.float64
else:
    bits = 32
    np_type = np.float32



@cuda.jit(device=True)
def gpu_euclidean_distance(X, Y):
    """Compute the euclidean distance between two vectors on GPU"""
    distance = 0.0
    for i in range(X.shape[0]):
        distance += (X[i] - Y[i]) ** 2
    return distance


@cuda.jit
def gpu_find_minimum_distance(mat, result_d, result_idx):
    """Compute distance between all rows of matrix X, keep track of minimum distance, and ignore distance that is not
       less than the minimum on GPU"""
    # Thread id in a 1D block
    tx = cuda.threadIdx.x
    # Block id in a 1D grid
    ty = cuda.blockIdx.x
    # Block width, i.e., number of threads per block
    bw = cuda.blockDim.x
    # Compute flattened index inside the array
    idx = tx + ty * bw
    if idx < mat.shape[0]:
        min_d = np.inf
        track_idx = cuda.local.array(2, dtype=np_type)
        for j in range(idx + 1, mat.shape[0]):
            d = gpu_euclidean_distance(mat[idx], mat[j])
            if d < min_d:
                min_d = d
                track_idx[0] = idx
                track_idx[1] = j

        result_d[idx] = min_d
        result_idx[idx, 0] = track_idx[0]
        result_idx[idx, 1] = track_idx[1]


def hierarchical_cluster(X):
    clusters = [[i] for i in range(X.shape[0])]
    # iteration_times = []
    with tqdm(total=len(clusters) - 1) as pbar:
        while len(clusters) > 1:
            # start_time = timeit.default_timer()
            mean_mat = [np.mean(X[clusters[i]], axis=0) for i in range(len(clusters))]
            mat = np.stack(mean_mat, axis=0)

            result_idx = np.zeros((mat.shape[0], 2)).astype(np.int32)
            result_d = np.zeros((mat.shape[0],)).astype(np.float32)
            threadsperblock = 32
            blockspergrid = (mat.shape[0] + (threadsperblock - 1)) // threadsperblock
            gpu_find_minimum_distance[blockspergrid, threadsperblock](mat, result_d, result_idx)
            min_distance = np.min(result_d)
            min_idx = np.argmin(result_d)
            track_idx = result_idx[min_idx]

            clusters[track_idx[0]].extend(clusters[track_idx[1]])
            clusters.pop(track_idx[1])

            # end_time = timeit.default_timer()  # End timer
            # time_taken = end_time - start_time
           # iteration_times.append(time_taken)
           #  print(len(clusters), time_taken)
            pbar.update(1)

    return np.array(clusters[0])