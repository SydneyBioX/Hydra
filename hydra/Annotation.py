#!/usr/bin/python

##############################################

# Manoj M Wagle (USydney, MIT CSAIL)

##############################################


import os, sys
import argparse
import logging
import glob
from tqdm import tqdm
from captum.attr import *
import scipy.special

import pandas as pd
import numpy as np
import pickle

# Create argument parser
parser = argparse.ArgumentParser("Annotation")

# Add arguments to the parser
parser.add_argument('--seed', type = int, default = 42, help ='seed')

# Input Data
parser.add_argument('--proccessedata_dir', metavar = 'DIR', default = os.getcwd(), help = 'Path to the directory containing processed single-cell data. Default: Current working directory')  
parser.add_argument('--balanceddata_dir', metavar = 'DIR', default = os.getcwd(), help = 'Path to the directory containing balanced single-cell data. Default: Current working directory')
parser.add_argument('--fs_dir', metavar = 'DIR', default = os.getcwd(), help = 'Path to the directory containing feature selection results. Default: Current working directory')
parser.add_argument('--modality', default='rna', choices=['rna', 'adt', 'atac'], help='Input data modality. Default: `rna`')

# For training
parser.add_argument('--batch_sizeanno', type = int, default = 64, help = 'batch size')
parser.add_argument('--epochsanno', type = int, default = 5, help = 'num of training epochs')
parser.add_argument('--lranno', type = float, default = 0.01, help = 'init learning rate')

# GPU specification    
parser.add_argument('--gpus', type = str, default = '0', help = 'Please specify the GPU to use')   

# Model
parser.add_argument('--num_classifiers', type = int, default = 25, help= 'Number of classifiers - Should be same as the number of models used for feature selection')

# Parse the command-line arguments
args = parser.parse_args()

if args.gpus:
    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpus


##############################################

# Import torch-related libraries after setting the CUDA_VISIBLE_DEVICES

import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from torch.utils.data import DataLoader, TensorDataset
from hydra.util import (setup_seed, real_label, NN_Classifier,
                  process_data_annotation_train, process_data_annotation_query)

device = torch.device("cuda" if torch.cuda.is_available()
                      else "mps" if getattr(torch, 'has_mps', False)
                      else "cpu")

print("Device to be used: ", device, "\n")

FloatTensor = torch.FloatTensor
LongTensor = torch.LongTensor

setup_seed(args.seed) ### set seed for reproducbility


##############################################

def main():

    dataset_folder_path = os.path.join(args.proccessedata_dir, "Input_Processed")

    if not os.path.exists(dataset_folder_path):
        logging.info(f"Error: The dataset folder '{dataset_folder_path}' does not exist. Please provide path to the directory containing `Input_Processed` directory!!!")
        sys.exit(1)

    dataset_folders = sorted(glob.glob(f"{args.proccessedata_dir}/Input_Processed"))

    for dataset_folder in dataset_folders:
        split_folders = sorted(glob.glob(f"{dataset_folder}/split_*"))

        for split_folder in split_folders:

            cwd = os.getcwd()
            feature_dir = f"{args.fs_dir}/Results/Feature_Selection/Hydra-{args.num_classifiers}/"

            args.rna = f"{split_folder}/rna_train.h5"
            if os.path.isfile(f"{split_folder}/adt_train.h5"):
                args.adt = f"{split_folder}/adt_train.h5"
            else:
                args.adt = "NULL" 
            if os.path.isfile(f"{split_folder}/atac_train.h5"):
                args.atac = f"{split_folder}/atac_train.h5"
            else:
                args.atac = "NULL"
            args.cty = f"{split_folder}/ct_train.csv"

            # Get a list of all balanced datasets in the split_folder
            balanced_data_files = sorted(glob.glob(f"{args.balanceddata_dir}/Balanced_Data-{args.num_classifiers}/train_data_bal_*.pt"))
            balanced_label_files = sorted(glob.glob(f"{args.balanceddata_dir}/Balanced_Data-{args.num_classifiers}/train_label_bal_*.pt"))

            model_save_dir = os.path.join(cwd, f"trained_model/Annotation_models-{args.num_classifiers}")
            os.makedirs(model_save_dir, exist_ok=True)  

            trainL = pd.read_csv(f"{split_folder}/ct_train.csv")

            save_path = os.path.join(model_save_dir, "ct_train.csv")

            trainL.to_csv(save_path, index=False)

            for idx, (data_file, label_file) in enumerate(zip(balanced_data_files, balanced_label_files)):

                # Train the classifer for the current balanced dataset
                logging.info(f"Training classifier: {idx+1}")

                args.balanced_data = data_file
                args.balanced_label = label_file
                
                (train_data, train_label, classify_dim) = process_data_annotation_train(args, feature_dir, split_folder)

                dataset = TensorDataset(train_data, train_label)
                dataloader = DataLoader(dataset, batch_size=args.batch_sizeanno, shuffle=True, drop_last = True)
                
                # Create NN classifier
                model = NN_Classifier(train_data.shape[1], classify_dim)

                if torch.cuda.device_count() > 1:
                    model = model.to(device)
                    model = nn.DataParallel(model)
                else:
                    model = model.to(device)

                criterion = nn.CrossEntropyLoss()
                optimizer = Adam(model.parameters(), args.lranno)

                # Training the model
                model.train()

                for epoch in tqdm(range(args.epochsanno)): 
                    epoch_loss = 0
                    correct_predictions = 0
                    total_predictions = 0

                    for batch in dataloader:
                        data, target = batch

                        data, target = data.to(device), target.to(device)

                        optimizer.zero_grad()

                        output = model(data)

                        loss = criterion(output, target)
                        epoch_loss += loss.item() 

                        loss.backward()

                        optimizer.step()

                        # Computing accuracy
                        predicted = torch.argmax(output, dim=1).to(device)
                        correct = (predicted == target).sum()  
                        correct_predictions += correct.item()
                        total_predictions += target.shape[0] 

                    epoch_loss /= len(dataloader)  # Taking average loss per epoch
                    accuracy = correct_predictions / total_predictions  # calculating accuracy

                    # print(f'Epoch {epoch+1}/{args.epochsanno}: Loss={epoch_loss}, Accuracy={accuracy}')

                # Save the model for this balanced dataset
                model_save_path = os.path.join(model_save_dir, f'FFNN_model_Hydra_{idx+1}.pt')
                final_model = model.module if isinstance(model, nn.DataParallel) else model
                torch.save(final_model.state_dict(), model_save_path)


            ##### Query dataset annotation #####

            logging.info("Annotating the query dataset...")
                    
            model_files = sorted(glob.glob(os.path.join(model_save_dir, "FFNN_model_Hydra_*.pt")))

            trainL = pd.read_csv(os.path.join(model_save_dir, "ct_train.csv"), index_col=0)

            known_cell_types = set(trainL['x'].values)

            if os.path.isfile(f"{split_folder}/adt_test.h5"):
                args.adt = f"{split_folder}/adt_test.h5"
            else:
                args.adt = "NULL" 
            if os.path.isfile(f"{split_folder}/atac_test.h5"):
                args.atac = f"{split_folder}/atac_test.h5"
            else:
                args.atac = "NULL"
            if os.path.isfile(f"{split_folder}/rna_test.h5"):
                args.rna = f"{split_folder}/rna_test.h5"
            else:
                args.rna = "NULL"

            (test_data, _, classify_dim) = process_data_annotation_query(args, feature_dir, split_folder)

            with open(os.path.join(split_folder, 'train_encoding.pkl'), 'rb') as f:
                train_encoding = pickle.load(f)
            
            transform_real_label, _ = real_label(args.cty, classify_dim)

            model_probabilities = []

            for model_file in model_files:
                model_load_path = model_file
                model = NN_Classifier(test_data.shape[1], len(known_cell_types))
                model.load_state_dict(torch.load(model_load_path, weights_only=False)) 
                model = model.to(device)
                model.eval()

                # Predict the labels for the test data
                with torch.no_grad():
                    output = model(test_data.to(device))
                    model_probabilities.append(output.cpu().numpy())  # For averaging, probability-based

            # Averaging the predicted probabilities from all models
            average_probabilities = np.mean(np.array(model_probabilities), axis=0)

            # Convert logits to probabilities
            average_probabilities = scipy.special.softmax(average_probabilities, axis=1)

            # Getting the class with the highest average probability
            test_pred = np.argmax(average_probabilities, axis=1)

            unique_cell_types = np.unique(transform_real_label)  

            unique_cell_types_dict = {}
            for cell_type in unique_cell_types:
                unique_cell_types_dict[cell_type] = train_encoding.get(cell_type, -1)
            
            # Reverse mapping from index to cell type name for predicted labels
            index_to_cell_type = {v: k for k, v in unique_cell_types_dict.items()}
            predicted_labels = [index_to_cell_type[idx] for idx in test_pred]

            # Creating a new DataFrame for the predicted labels similar to the ground truth format
            predicted_labels_df = pd.DataFrame(predicted_labels, columns=['x'])
            predicted_labels_df.index += 1  # Adjust index to start from 1 to match the example format
            predicted_labels_df.index.name = ""  # Match the format which has no name for index column

            # Saving annotation results
            cwd = os.getcwd()  # Get the current working directory
            Hydra_folder = os.path.join(cwd, f"Results/Annotation/Hydra-{args.num_classifiers}")
            if not os.path.exists(Hydra_folder):
                os.makedirs(Hydra_folder)

            # Save the predicted labels DataFrame as a CSV file
            predicted_labels_path = os.path.join(Hydra_folder, f"cell_type_predicted_Hydra-{args.num_classifiers}.csv")
            predicted_labels_df.to_csv(predicted_labels_path)

    return


##############################################

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format = Log_Format, level = logging.INFO)

    # Call the main function
    main()
