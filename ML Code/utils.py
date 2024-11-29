import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import torch
from LSTM import EstimationLSTM

class Utils:


    @staticmethod
    def plot_outputs(outputs, labels, min_val, max_val):
        outputs_denorm = Utils.denormalization(outputs, min_val, max_val)
        labels_denorm = Utils.denormalization(labels, min_val, max_val)
        print("Denormalized Outputs:")
        print(outputs_denorm[0:10, :])
        print("Denormalized Labels:")
        print(labels_denorm[0:10, :])
        plt.plot(outputs_denorm, 'r', labels_denorm, 'b')
        # plt.show()
    @staticmethod
    def make_dir(path):
        os.makedirs(path, exist_ok=True)

    @staticmethod
    def save_log(log_dir,loss_log):
        logs =np.array(loss_log)
        log_df = pd.DataFrame(logs,columns=["Test",'Validation'])
        log_df.to_csv(log_dir+'/log.csv')




class ModelManager:
    def __init__(self, device):
        self.device = device

    def model_setting(self, input_size, hidden_size, output_size, seq_len, num_layers):
        self.model = EstimationLSTM(input_size, hidden_size, output_size, seq_len, num_layers).to(self.device)
        return self.model

    def save_checkpoint(self, model, path, onnx_path, batch_size, seq_len, input_size):
        torch.save(model, path)
        dummy_input = torch.randn(batch_size, seq_len, input_size).to(device=self.device)
        torch.onnx.export(model, dummy_input, onnx_path, export_params=True, opset_version=10,
                          input_names=['Modelinput'], output_names=['Modeloutput'])

    def load_checkpoint(self, path):
        self.model = torch.load(path, map_location=self.device)
        return self.model
    