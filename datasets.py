from torchvision.datasets import MNIST
import matplotlib.pyplot as plt
import time
import torch
import torch.nn as nn
import torch.optim as optim
from torch.distributions import MultivariateNormal
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
from scipy.signal import medfilt
import torch.nn.functional as F
from sklearn.metrics import mean_squared_error
from torch.utils.data import random_split

def transform_img2pc(img,thresh=120):
    # transform image to point cloud
    # TODO: replace constant threshold to probablistsic value
    img_array = np.asarray(img)
    indices = np.argwhere(img_array > thresh )
    return indices.astype(np.float32)

class MNIST2D(Dataset):
    """MNIST pointcloud dataset."""

    def __init__(self, dataset, num_points):
        self.dataset = dataset
        self.number_of_points = num_points

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, idx):

        img,label = self.dataset[idx]
        pc = transform_img2pc(img)
        
        sampling_indices = np.random.choice(pc.shape[0], min(len(pc),self.number_of_points),replace=False)
        pc = pc[sampling_indices, :]


        pc = pc.astype(np.float32) 
        pc = torch.tensor(pc)

        return pc, label