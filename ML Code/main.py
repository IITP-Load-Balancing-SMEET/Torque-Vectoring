import os
from datetime import datetime
from dataloader import DataHandler
from utils import Utils, ModelManager
from trainer import Trainer
import torch
import numpy as np

if __name__ == "__main__":
    try:
        batch_size = 8
        seq_len = 500
        hidden_sizes = [100]
        epochs = 200
        input_size = 15
        num_layers = 3
        pretrained = False

        x_min_value = np.array([-0.6, 0, 0, 0, 0, 0, 0, 0, 0, -30, -30,-20,-1, -2, 0])
        x_max_value = np.array([0.6, 30, 30, 30, 30,1700 ,1700 ,1700 ,1700, 30,30,20,1, 2, 35])
        y_min_value = np.array([-3000,-3000,-3000,-3000])
        y_max_value = np.array([3000,3000,3000,3000])

        today = datetime.today()
        data_path = os.path.join('Data', 'Test')
        log_dir = os.path.join('Data/ML', str(today.date()), 'log')
        pt_dir = os.path.join('Data/ML', str(today.date()), "pt")
        onnx_dir = os.path.join('Data/ML', str(today.date()), "onnx")
        pt_path = 'Data/ML/pretrained_model/H_40_FC_1.pt'

        for hidden_size in hidden_sizes:
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            model_manager = ModelManager(device)
            if pretrained:
                model = model_manager.load_checkpoint(pt_path)
            else:
                model = model_manager.model_setting(input_size, hidden_size, output_size=4, seq_len=seq_len, num_layers=num_layers)
            trainer = Trainer(model, device, data_path,batch_size,input_size,seq_len,x_max_value,x_min_value,y_max_value,y_min_value)
            Utils.make_dir(log_dir)
            Utils.make_dir(pt_dir)
            Utils.make_dir(onnx_dir)
            loss_log = trainer.train(epochs)
            model_manager.save_checkpoint(model, os.path.join(pt_dir, f'model_epoch_{epochs}.pt'),
                                          os.path.join(onnx_dir, f'model_epoch_{epochs}.onnx'),
                                          batch_size, seq_len, input_size)
            # Save logs
            Utils.save_log(log_dir, loss_log)
    except KeyboardInterrupt:
        print("Canceled by user...")