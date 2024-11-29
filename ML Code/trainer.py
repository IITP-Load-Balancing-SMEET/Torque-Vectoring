from tqdm import tqdm
import torch
import torch.nn as nn
import torch.optim as optim
from utils import Utils
from dataloader import DataHandler
import os
class Trainer:
    def __init__(self, model, device, data_path, batch_size, input_size, seq_len, x_max_value,x_min_value,
                 y_max_value,y_min_value,loss = nn.MSELoss()
                ):
        self.model = model
        self.device = device
        self.data_path = data_path

        self.batch_size = batch_size
        self.input_size = input_size
        self.seq_len = seq_len

        self.x_max_value = x_max_value
        self.y_max_value = y_max_value

        self.x_min_value = x_min_value
        self.y_min_value = y_min_value

        self.criterion = loss
        self.optimizer = optim.Adam(self.model.parameters(), lr=1e-4)
        self.loss_log = []
        self.data_handler = DataHandler(self.device, x_max_value,x_min_value,
                 y_max_value,y_min_value,)

  
    def train_epoch(self,cu_train_loss):
        total_train_loss = 0
        self.model.train()

        for inputs, labels in self.train_loader:
            self.model.reset_hidden_state()
            self.optimizer.zero_grad()
            outputs = self.model(inputs)
                     
            # Ensure the outputs and labels have the same shape
            labels = labels.squeeze(dim=1)
            train_loss = self.criterion(outputs, labels)
            train_loss.backward()
            self.optimizer.step()
            total_train_loss += train_loss.item()
        avg_train_loss = total_train_loss / len(self.train_loader)
        return (avg_train_loss+cu_train_loss)/2

    def validate(self,cu_val_loss):
        self.model.eval()
        val_loss = 0
        with torch.no_grad():
            for inputs, labels in self.val_loader:
                outputs = self.model(inputs)
                labels = labels.squeeze(dim=1)
                val_loss += self.criterion(outputs, labels).item()
        avg_val_loss = val_loss / len(self.val_loader)


        return (avg_val_loss+cu_val_loss)/2


    def train(self, epochs):
        for epoch in tqdm(range(epochs)):
            if epoch in [50, 100, 200]:
                for param_group in self.optimizer.param_groups:
                    param_group['lr'] *= 0.1
            
        # import data part
            data_list = self.data_handler.make_data_list(self.data_path)
            train_loss = 0
            val_loss = 0
            for data_name in data_list:
                print("Now Proccesing :" + data_name)
                current_data = os.path.join(self.data_path,data_name)
                self.train_loader, self.val_loader, _ = self.data_handler.data_loader(current_data,self.batch_size,self.seq_len)
                train_loss = self.train_epoch(train_loss)
                val_loss = self.validate(val_loss)
            self.loss_log.append([train_loss, val_loss])
            print(f'Epoch {epoch + 1}, Loss: {train_loss:.8f}, Validation Loss: {val_loss:.8f}')
        return self.loss_log


        # # check_the data
        # test_inputs, test_labels = next(iter(self.val_loader))
        # outputs = self.model(test_inputs)
        # outputs = outputs.cpu().detach().numpy()
        # test_labels = test_labels.cpu().detach().numpy()
        # outputs = self.data_handler.denormalization(outputs,self.y_min_value,self.y_max_value)
        # test_labels = self.data_handler.denormalization(test_labels,self.y_min_value,self.y_max_value)
        # print("Estimation")
        # print(outputs)
        # print("Ground_Truth")
        # print(test_labels)