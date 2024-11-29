import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, random_split
from dataset import CustomDataset
from utils import Utils

class DataHandler:
    def __init__(self, device, x_max_value,x_min_value,
                y_max_value,y_min_value,):
        self.device = device
        self.x_max_value = x_max_value
        self.y_max_value = y_max_value

        self.x_min_value = x_min_value
        self.y_min_value = y_min_value

    def make_data_list(self, data_path):
        file_list = os.listdir(data_path)

        return file_list

    def import_data(self, seq_len, data_path):
        data_X = []
        data_Y = []
        data = pd.read_csv(data_path)
        y = data.iloc[:,-4:].values
        data = data.iloc[:,:-4].values
        for i in range(0, len(data) - seq_len):
            _x = data[i:i + seq_len, :]
            _y = np.array([y[i + seq_len - 1]])
            data_X.append(_x)
            data_Y.append(_y)
        return data_X, data_Y

    def data_loader(self, data_path, batch_size=32, seq_len=200, split=0.8):
        self.seq_len = seq_len
        self.batch_size = batch_size
        data_X, data_Y = self.import_data(seq_len=self.seq_len, data_path=data_path)
        self.np_x = np.array(data_X)
        self.np_y = np.array(data_Y)
        norm_x = self.normalization(self.np_x, self.x_min_value, self.x_max_value)
        norm_y = self.normalization(self.np_y, self.y_min_value, self.y_max_value)
        data_X = torch.FloatTensor(norm_x).to(self.device)
        data_Y = torch.FloatTensor(norm_y).to(self.device)
        dataset = CustomDataset(data_X, data_Y)
        train_dataset, val_dataset = random_split(dataset, [int(len(data_X) * split), len(data_X) - int(len(data_X) * split)])
        self.train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        self.val_loader = DataLoader(val_dataset, batch_size=batch_size)
        sample_batch = next(iter(self.train_loader))
        self.input_size = sample_batch[0].shape[2]
        return self.train_loader, self.val_loader, self.input_size
    
    def check_data(self,data):
        denorm_y = Utils.denormalization(data, self.y_min_value,self.y_max_value)
        return denorm_y
    
    def normalization(self,data, min_val, max_val):
        data_shape = data.shape
        data = np.reshape(data, (-1, data.shape[-1]))
        min_value = np.zeros(shape=np.size(min_val))
        max_value = np.ones(shape=np.size(max_val))
        scaled_data = min_value + (max_value - min_value) * (data - min_val) / (max_val - min_val)
        scaled_data = np.clip(scaled_data, 0, 1)
        scaled_data = np.reshape(scaled_data, data_shape)
        return scaled_data

    def denormalization(self,scaled_data, min_val, max_val):
        min_val = np.array(min_val).reshape(1, 1, -1)  # Reshape to [1, 1, feature_dim]
        max_val = np.array(max_val).reshape(1, 1, -1)  # Reshape to [1, 1, feature_dim]
        denormalized_data = min_val + (max_val - min_val) * scaled_data
        return denormalized_data