#!/usr/bin/python

##############################################

# Manoj M Wagle (USydney; CMRI)

##############################################


import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
import scanpy as sc
import random, os, re
from torch.autograd import Variable
import h5py
import scipy
import pickle

device = torch.device("cuda" if torch.cuda.is_available()
                      else "mps" if getattr(torch, 'has_mps', False)
                      else "cpu")

FloatTensor = torch.FloatTensor
LongTensor = torch.LongTensor


##################################################

def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     torch.backends.cudnn.benchmark = False
     os.environ['PYTHONHASHSEED']=str(seed)


##################################################

def real_label(label_path, classify_dim):
    output_v = []
    output_dict = {}
    label = pd.read_csv(label_path, header=None, index_col=False)
    label_real = label.iloc[1:(label.shape[0]), 1]
    (label_num, _) = read_fs_label(label_path)
    
    # Get unique labels in sorted order
    unique_labels = sorted(np.unique(label_num.cpu()))

    for i in unique_labels:    
        idx = np.array(torch.where(label_num==i)[0].cpu()).astype('int32')
        if len(idx) > 0:
            temp = label_real[idx[0] + 1]
            output_v.append(temp)
            count = sum(label_real == temp)  # Count the occurrences of temp in label_real
            output_dict[i] = (temp, count)
    return output_v, output_dict


##################################################

def Index2Label(label_path,classify_dim):
    label = pd.read_csv(label_path,header=None,index_col=False)  
    label_real = label.iloc[1:(label.shape[0]),1]
    (label_num, _) = read_fs_label(label_path)
    index_to_label = {}
    for i in range(classify_dim):    
        temp = label_real[np.array(torch.where(label_num==i)[0][0].cpu()).astype('int32')+1]
        index_to_label[i] = temp
    return index_to_label


##################################################

def read_fs_label(label_path):
    label_fs = pd.read_csv(label_path, header=None, index_col=False)
    label_fs_codes = pd.Categorical(label_fs.iloc[1:(label_fs.shape[0]),1]).codes
    label_fs_names = label_fs.iloc[1:, 1].tolist()
    # for x,y in zip(label_fs_names[:15], label_fs_codes[:15]):
    #     print(x, " : ", y)
    # sys.exit()
    label_fs_codes = np.array(label_fs_codes[:]).astype('int32')
    label_fs_codes = torch.from_numpy(label_fs_codes)
    label_fs_codes = label_fs_codes.type(LongTensor)
    return label_fs_codes, label_fs_names


##################################################

def read_fs_label_test(label_path, encoding = None):
    label_fs = pd.read_csv(label_path, header=None, index_col=False, skiprows=1)
    label_fs_names = label_fs[1].tolist()
    label_fs['encoded'] = label_fs[1].map(encoding)
    label_fs_codes = torch.from_numpy(label_fs['encoded'].values).type(LongTensor)
    # for x,y in zip(label_fs_names[:15], label_fs_codes[:15]):
    #     print(x, " : ", y)
    # sys.exit()
    return label_fs_codes, label_fs_names


##################################################

class MyDataset(Dataset):
    def __init__(self, data, label):
        self.data = data
        self.label = label
        self.labels = label

    def __getitem__(self, index): 
        img, target = self.data[index,:], self.label[index]
        sample = {'data': img, 'label': target}
        return sample

    def __len__(self):
        return len(self.data)
        

##################################################

class ToTensor(object):
    def __call__(self, sample):
        data,label = sample['data'], sample['label']
        data = data.transpose((1, 0))
        return {'data': torch.from_numpy(data),
                'label': torch.from_numpy(label),
               }
    

##################################################
               
class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self, name, fmt=':f'):
        self.name = name
        self.fmt = fmt
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

    def __str__(self):
        fmtstr = '{name} {val' + self.fmt + '} ({avg' + self.fmt + '})'
        return fmtstr.format(**self.__dict__)


##################################################

def save_checkpoint(state, filename):
    save_dir = os.path.dirname(filename)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    torch.save(state, filename)
    

##################################################

class CrossEntropyLabelSmooth(nn.Module):

  def __init__(self, num_classes=17, epsilon=0.1):
    super(CrossEntropyLabelSmooth, self).__init__()
    self.num_classes = num_classes
    self.epsilon = epsilon
    self.logsoftmax = nn.LogSoftmax(dim=1)

  def forward(self, inputs, targets):
    log_probs = self.logsoftmax(inputs)
    targets = torch.zeros_like(log_probs).scatter_(1, targets.unsqueeze(1), 1)
    targets = (1 - self.epsilon) * targets + self.epsilon / self.num_classes
    loss = (-targets * log_probs).mean(0).sum()
    return loss


##################################################

def read_h5_data(data_path):
    data = h5py.File(data_path,"r")
    h5_data = data['matrix/data']
    sparse_data = scipy.sparse.csr_matrix(np.array(h5_data).transpose())
    data_fs = torch.from_numpy(np.array(sparse_data.todense()))
    data_fs = Variable(data_fs.type(FloatTensor))
    return data_fs


##################################################

def get_encodings(model, dl):
    model.eval()
    encodings = []
    ori_data = []
    label = [] 
    with torch.no_grad():
        for i, batch_sample in enumerate(dl):
                x = batch_sample['data']
                x = Variable(x.type(FloatTensor))
                x = torch.reshape(x,(x.size(0),-1))
                rna_valid_label = batch_sample['label']
                rna_valid_label = Variable(rna_valid_label.type(LongTensor))
            
                x_prime = model.encoder(x.to(device))
                encodings.append(x_prime)
                ori_data.append(x)
                label.append(rna_valid_label)
                
                
    return torch.cat(encodings, dim=0),torch.cat(ori_data,dim=0),torch.cat(label,dim=0)


##################################################

def get_decodings(model, dl):
    model.eval()
    decodings = []
    ori_data = []
    with torch.no_grad():
        for i, batch_sample in enumerate(dl):
                x = batch_sample['data']
                x = Variable(x.type(FloatTensor))
                x = torch.reshape(x,(x.size(0),-1))
                rna_valid_label = batch_sample['label']
                rna_valid_label = Variable(rna_valid_label.type(LongTensor))
                x_prime, x_cty,mu,var = model(x.to(device))
                decodings.append(x_prime)
                ori_data.append(x)                
    return torch.cat(decodings, dim=0),torch.cat(ori_data,dim=0)


##################################################

def get_simulated_data_random_generation(model, dl):
    model.eval()
    decodings = []
    ori_data = []
    label = []
    for i, batch_sample in enumerate(dl):
            x = batch_sample['data']
            x = Variable(x.type(FloatTensor))
            x = torch.reshape(x,(x.size(0),-1))
            rna_valid_label = batch_sample['label']
            rna_valid_label = Variable(rna_valid_label.type(LongTensor))
            
            x_prime, x_cty,mu,var = model(x.to(device))
            decodings.append(x_prime)
            ori_data.append(x)
            label.append(rna_valid_label)
                
    return torch.cat(decodings, dim=0),torch.cat(label,dim=0)


##################################################

def get_simulated_data(model, dl):
    model.eval()
    decodings = []
    ori_data = []
    label = []
    with torch.no_grad():
        for i, batch_sample in enumerate(dl):
                x = batch_sample['data']
                x = Variable(x.type(FloatTensor))
                x = torch.reshape(x,(x.size(0),-1))
                rna_valid_label = batch_sample['label']
                rna_valid_label = Variable(rna_valid_label.type(LongTensor))
            
                x_prime, x_cty,mu,var = model(x.to(device))
                decodings.append(x_prime)
                ori_data.append(x)
                label.append(rna_valid_label)
                
    return torch.cat(decodings, dim=0),torch.cat(label,dim=0),torch.cat(ori_data,dim=0)
    

##################################################

def KL_loss(mu, logvar):
    KLD = -0.5 * torch.mean(1 + logvar - mu**2 -  logvar.exp())
    return  KLD


##################################################

def get_vae_simulated_data_from_sampling(model, dl):
    model.eval()
    latent = []
    ori_data = []
    label = []
    decodings = []
    with torch.no_grad():
        for i, batch_sample in enumerate(dl):
                x = batch_sample['data']
                x = Variable(x.type(FloatTensor))
                x = torch.reshape(x,(x.size(0),-1))
                
                ori_data.append(x)
                rna_valid_label = batch_sample['label']
                rna_valid_label = Variable(rna_valid_label.type(LongTensor))
                y, x_cty,mu,var = model(x.to(device))
                decodings.append(y)
                label.append(rna_valid_label)
                
    return torch.cat(decodings, dim=0),torch.cat(label,dim=0),torch.cat(ori_data,dim=0)


##############################################

def load_and_preprocess_data(args, setting=""):
    if args.adt != "NULL" and args.atac != "NULL":

        mode = "TEAseq"

        if setting == 'train':
            train_rna_data_path = args.rna
            train_adt_data_path = args.adt
            train_atac_data_path = args.atac
            train_label_path = args.cty

            train_rna_data = read_h5_data(train_rna_data_path)
            train_adt_data = read_h5_data(train_adt_data_path)
            train_atac_data = read_h5_data(train_atac_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(train_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = train_rna_data.shape[1]
            nfeatures_adt = train_adt_data.shape[1]
            nfeatures_atac = train_atac_data.shape[1]
            
            # Total features count 
            feature_num = nfeatures_rna + nfeatures_adt + nfeatures_atac

            # Prepare training data
            train_data = torch.cat((train_rna_data, train_adt_data, train_atac_data), 1) 
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_rna_data_path = args.rna
            test_adt_data_path = args.adt
            test_atac_data_path = args.atac
            test_label_path = args.cty

            test_rna_data = read_h5_data(test_rna_data_path)
            test_adt_data = read_h5_data(test_adt_data_path)
            test_atac_data = read_h5_data(test_atac_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(test_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = test_rna_data.shape[1]
            nfeatures_adt = test_adt_data.shape[1]
            nfeatures_atac = test_atac_data.shape[1]

            # Total features count 
            feature_num = nfeatures_rna + nfeatures_adt + nfeatures_atac

            # Prepare test data
            test_data = torch.cat((test_rna_data, test_adt_data, test_atac_data), 1) 
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)


    if args.adt != "NULL" and args.atac == "NULL":

        mode = "CITEseq"

        if setting == 'train':
            train_rna_data_path = args.rna
            train_adt_data_path = args.adt
            train_label_path = args.cty

            train_rna_data = read_h5_data(train_rna_data_path)
            train_adt_data = read_h5_data(train_adt_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(train_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = train_rna_data.shape[1]
            nfeatures_adt = train_adt_data.shape[1]
            
            # Total features count 
            feature_num = nfeatures_rna + nfeatures_adt

            # Prepare training data
            train_data = torch.cat((train_rna_data, train_adt_data), 1) 
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_rna_data_path = args.rna
            test_adt_data_path = args.adt
            test_label_path = args.cty

            test_rna_data = read_h5_data(test_rna_data_path)
            test_adt_data = read_h5_data(test_adt_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(test_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = test_rna_data.shape[1]
            nfeatures_adt = test_adt_data.shape[1]

            # Total features count 
            feature_num = nfeatures_rna + nfeatures_adt

            # Prepare test data
            test_data = torch.cat((test_rna_data, test_adt_data), 1) 
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)


    if args.adt == "NULL" and args.atac != "NULL":

        mode = "SHAREseq"

        if setting == 'train':
            train_rna_data_path = args.rna
            train_atac_data_path = args.atac
            train_label_path = args.cty

            train_rna_data = read_h5_data(train_rna_data_path)
            train_atac_data = read_h5_data(train_atac_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(train_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = train_rna_data.shape[1]
            nfeatures_atac = train_atac_data.shape[1]

            # Total features count 
            feature_num = nfeatures_rna + nfeatures_atac

            # Prepare training data
            train_data = torch.cat((train_rna_data, train_atac_data), 1) 
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_rna_data_path = args.rna
            test_atac_data_path = args.atac
            test_label_path = args.cty

            test_rna_data = read_h5_data(test_rna_data_path)
            test_atac_data = read_h5_data(test_atac_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(test_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = test_rna_data.shape[1]
            nfeatures_atac = test_atac_data.shape[1]

            # Total features count 
            feature_num = nfeatures_rna + nfeatures_atac

            # Prepare test data
            test_data = torch.cat((test_rna_data, test_atac_data), 1) 
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)


    if args.adt == "NULL" and args.atac == "NULL":

        mode = "scRNA-Seq"

        if setting == 'train':
            train_rna_data_path = args.rna
            train_label_path = args.cty

            train_rna_data = read_h5_data(train_rna_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(train_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = train_rna_data.shape[1]

            # Total features count 
            feature_num = nfeatures_rna

            # Prepare training data
            train_data = train_rna_data
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_rna_data_path = args.rna
            test_label_path = args.cty

            test_rna_data = read_h5_data(test_rna_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()

            _, label_to_name_mapping = real_label(test_label_path, classify_dim)
            
            # Get the number of features
            nfeatures_rna = test_rna_data.shape[1]
            
            # Total features count 
            feature_num = nfeatures_rna

            # Prepare test data
            test_data = test_rna_data
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)
    
    if args.rna == "NULL" and args.atac == "NULL":

        mode = "scADT-Seq"

        if setting == 'train':
            train_adt_data_path = args.adt
            train_label_path = args.cty
            train_adt_data = read_h5_data(train_adt_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()
            _, label_to_name_mapping = real_label(train_label_path, classify_dim)

            # Get the number of features
            nfeatures_adt = train_adt_data.shape[1]
            
            # Total features count 
            feature_num = nfeatures_adt
            
            # Prepare training data
            train_data = train_adt_data
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_adt_data_path = args.adt
            test_label_path = args.cty
            test_adt_data = read_h5_data(test_adt_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()
            _, label_to_name_mapping = real_label(test_label_path, classify_dim)

            # Get the number of features
            nfeatures_adt = test_adt_data.shape[1]

            # Total features count 
            feature_num = nfeatures_adt
            
            # Prepare test data
            test_data = test_adt_data
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)
    

    if args.rna == "NULL" and args.adt == "NULL":

        mode = "scATAC-Seq"
    
        if setting == 'train':
            train_atac_data_path = args.atac
            train_label_path = args.cty
            train_atac_data = read_h5_data(train_atac_data_path)
            (train_label, _) = read_fs_label(train_label_path)
            classify_dim = (max(train_label)+1).cpu().numpy()
            _, label_to_name_mapping = real_label(train_label_path, classify_dim)

            # Get the number of features
            nfeatures_atac = train_atac_data.shape[1]
            
            # Total features count 
            feature_num = nfeatures_atac
            
            # Prepare training data
            train_data = train_atac_data
            train_transformed_dataset = MyDataset(train_data, train_label)
            train_dl = DataLoader(train_transformed_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)
        else:
            test_atac_data_path = args.atac
            test_label_path = args.cty
            test_atac_data = read_h5_data(test_atac_data_path)
            (test_label, _) = read_fs_label(test_label_path)
            classify_dim = (max(test_label)+1).cpu().numpy()
            _, label_to_name_mapping = real_label(test_label_path, classify_dim)

            # Get the number of features
            nfeatures_atac = test_atac_data.shape[1]

            # Total features count 
            feature_num = nfeatures_atac
            
            # Prepare test data
            test_data = test_atac_data
            test_transformed_dataset = MyDataset(test_data, test_label)
            test_dl = DataLoader(test_transformed_dataset, batch_size = args.batch_size, shuffle = False, num_workers = 0, drop_last = False)


    if mode == "TEAseq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_rna, nfeatures_adt, nfeatures_atac, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_rna, nfeatures_adt, nfeatures_atac, feature_num, label_to_name_mapping

    if mode == "CITEseq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_rna, nfeatures_adt, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_rna, nfeatures_adt, feature_num, label_to_name_mapping

    if mode == "SHAREseq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_rna, nfeatures_atac, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_rna, nfeatures_atac, feature_num, label_to_name_mapping

    if mode == "scRNA-Seq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_rna, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_rna, feature_num, label_to_name_mapping

    if mode == "scADT-Seq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_adt, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_adt, feature_num, label_to_name_mapping

    if mode == "scATAC-Seq":
        if setting == 'train':
            return train_data, train_dl, train_label, mode, classify_dim, nfeatures_atac, feature_num, label_to_name_mapping
        else:
            return test_data, test_dl, test_label, mode, classify_dim, nfeatures_atac, feature_num, label_to_name_mapping


##################################################

def process_data_annotation_train(args, feature_dir, split_folder):

    label_path = args.cty

    # Load data
    train_data = torch.load(args.balanced_data)
    train_label = torch.load(args.balanced_label)

    (label_fs_codes, label_fs_names) = read_fs_label(label_path)
    # train_encoding dictionary
    train_encoding = {}
    for name, code in zip(label_fs_names, label_fs_codes):
        if name not in train_encoding:
            train_encoding[name] = code.cpu().numpy().item()

    # Add the 'Unknown' key to train_encoding with the next available integer code
    train_encoding['Unknown'] = max(train_encoding.values()) + 1
    with open(os.path.join(split_folder, 'train_encoding.pkl'), 'wb') as f:
        pickle.dump(train_encoding, f)

    classify_dim = (max(train_label)+1).cpu().numpy()

    ### Hybrid Feature Selection Strategy
    # Load the original labels
    df_original = pd.read_csv(label_path)
    cell_type_counts = df_original['x'].value_counts()

    # Determine the threshold for rare and major classes
    median_count = cell_type_counts.median()
    rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

    top_features = []

    # Process each feature file
    for filename in os.listdir(feature_dir):
        if filename.endswith('.csv'):
            filepath = os.path.join(feature_dir, filename)
            
            # Read the CSV file
            df = pd.read_csv(filepath)
            
            # Get cell type from the filename
            cell_type = re.search('fs\.(.*?)_', filename).group(1)
            
            # number of top features to extract
            if cell_type in rare_cell_types:
                num_features = 100  # number of top features for minor cell types
            else:
                num_features = 100  # number of top features for major cell types
            
            # Get the top features
            top_features_file = df['Feature Name'].head(num_features).tolist()
            
            # Add these features to the list
            top_features.extend(top_features_file)

    # Get the unique set of top features across all files
    top_features = list(set(top_features))

    # Extract indices from the features and sort them
    top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

    # Filter data to only include top features
    train_data = train_data[:, top_feature_indices]

    return train_data, train_label, classify_dim


##############################################

def process_data_annotation_query(args, feature_dir, split_folder):
    if args.adt != "NULL" and args.atac != "NULL":
        mode = "TEAseq"

        rna_data_path = args.rna
        adt_data_path = args.adt
        atac_data_path = args.atac

        train_label = torch.load(args.balanced_label)

        # Load data
        rna_data = read_h5_data(rna_data_path)
        adt_data = read_h5_data(adt_data_path)
        atac_data = read_h5_data(atac_data_path)

        test_data = torch.cat((rna_data, adt_data, atac_data), 1)

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []

        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to extract
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for minor cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                top_features.extend(top_features_file)

        # Get the unique set of top features across all files
        top_features = list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]

    if args.adt != "NULL" and args.atac == "NULL":
        mode = "CITEseq"

        rna_data_path = args.rna
        adt_data_path = args.adt

        train_label = torch.load(args.balanced_label)

        # Load data
        rna_data = read_h5_data(rna_data_path)
        adt_data = read_h5_data(adt_data_path)

        test_data = torch.cat((rna_data, adt_data), 1)    

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []

        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to extract
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for minor cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                top_features.extend(top_features_file)

        # Get the unique set of top features across all files
        top_features = list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]

    if args.adt == "NULL" and args.atac != "NULL":
        mode = "SHAREseq"

        rna_data_path = args.rna
        atac_data_path = args.atac

        train_label = torch.load(args.balanced_label)

        # Load data
        rna_data = read_h5_data(rna_data_path)
        atac_data = read_h5_data(atac_data_path)

        test_data = torch.cat((rna_data, atac_data), 1)    

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []

        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to extract
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for minor cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                top_features.extend(top_features_file)

        # Get the unique set of top features across all files
        top_features = list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]

    if args.adt == "NULL" and args.atac == "NULL":

        mode = "scRNA-Seq"

        rna_data_path = args.rna

        train_label = torch.load(args.balanced_label)
        
        test_data = read_h5_data(rna_data_path)

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []
        
        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to select
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for rare cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                for feature in top_features_file:
                    if feature not in top_features:
                        top_features.append(feature)

        # Get the unique set of top features across all files
        top_features= list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]
    
    if args.rna == "NULL" and args.atac == "NULL":

        mode = "scADT-Seq"

        adt_data_path = args.adt

        train_label = torch.load(args.balanced_label)
    
        test_data = read_h5_data(adt_data_path)

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []
        
        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to select
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for rare cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                for feature in top_features_file:
                    if feature not in top_features:
                        top_features.append(feature)

        # Get the unique set of top features across all files
        top_features= list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]

    if args.rna == "NULL" and args.adt == "NULL":

        mode = "scATAC-Seq"

        atac_data_path = args.atac

        train_label = torch.load(args.balanced_label)
    
        test_data = read_h5_data(atac_data_path)

        classify_dim = (max(train_label)+1).cpu().numpy()

        ### Hybrid Feature Selection Strategy
        # Load the original labels
        df_original = pd.read_csv(args.cty)
        cell_type_counts = df_original['x'].value_counts()

        # Determine the threshold for rare and major classes
        median_count = cell_type_counts.median()
        rare_cell_types = cell_type_counts[cell_type_counts < median_count].index.tolist()

        top_features = []
        
        # Process each feature file
        for filename in os.listdir(feature_dir):
            if filename.endswith('.csv'):
                filepath = os.path.join(feature_dir, filename)
                
                # Read the CSV file
                df = pd.read_csv(filepath)
                
                # Get cell type from the filename
                cell_type = re.search('fs\.(.*?)_', filename).group(1)
                
                # number of top features to select
                if cell_type in rare_cell_types:
                    num_features = 100  # number of top features for rare cell types
                else:
                    num_features = 100  # number of top features for major cell types
                
                # Get the top features
                top_features_file = df['Feature Name'].head(num_features).tolist()
                
                # Add these features to the list
                for feature in top_features_file:
                    if feature not in top_features:
                        top_features.append(feature)

        # Get the unique set of top features across all files
        top_features= list(set(top_features))

        # Extract indices from the features and sort them
        top_feature_indices = sorted([int(feature.split("_")[-1]) for feature in top_features])

        # Filter data to only include top features
        test_data = test_data[:, top_feature_indices]


    return test_data, train_label, classify_dim


##############################################

class Mish(nn.Module):
    def forward(self, x):
        return x * torch.tanh(F.softplus(x))


##############################################

class NN_Classifier(nn.Module):
    def __init__(self, feature_num, num_classes):
        super(NN_Classifier, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(feature_num, 128), 
            Mish(),                        
            nn.Linear(128, num_classes)
        )
    def forward(self, x):
        return self.net(x)


##############################################

def accuracy(output, target, topk=(1,)):
    """Computes the accuracy over the k top predictions for the specified values of k"""
    with torch.no_grad():
        target = target.to(device)

        maxk = max(topk)
        batch_size = target.size(0)

        _, pred = output.topk(maxk, 1, True, True)
        pred = pred.t()
        correct = pred.eq(target.view(1, -1).expand_as(pred))

        res = []
        for k in topk:
            correct_k = correct[:k].view(-1).float().sum(0, keepdim=True)
            res.append(correct_k.mul_(100.0 / batch_size))
            
        return res


##################################################

def clip_counts(new_data, new_label, median_index_length):
    unique_labels = torch.unique(new_label)
    
    indices_to_keep = []
    
    for label in unique_labels:
        indices_of_samples = torch.nonzero(new_label == label, as_tuple=True)[0]
        
        if len(indices_of_samples) > median_index_length:
            # Randomly select samples to keep
            indices_to_keep += indices_of_samples[torch.randperm(len(indices_of_samples))[:median_index_length]].tolist()
        else:
            indices_to_keep += indices_of_samples.tolist()
            
    # Select desired samples from new_data and new_label
    new_data = new_data[indices_to_keep]
    new_label = new_label[indices_to_keep]
    
    return new_data, new_label


##################################################

def perform_data_augmentation(label_to_name_mapping, train_num, classify_dim, train_label, train_data, model, args):

    # Create an empty list to store tuples of class index and the number of samples in that class
    stage1_list = []

    new_label_names = []

    # Iterate through each class index and append a tuple to the list
    for i in np.arange(0, classify_dim):
        stage1_list.append([i, train_num[i]])
    stage1_df = pd.DataFrame(stage1_list)

    # Create a copy of stage1_df for adjustments
    adjusted_df = stage1_df.copy()

    # Calculate the median number of samples for classes
    sorted_train_num = np.sort(train_num)

    if classify_dim % 2 == 0:
        mid_idx = int(classify_dim / 2)
        middle_1 = sorted_train_num[mid_idx - 1]
        middle_2 = sorted_train_num[mid_idx]
        
        # Calculate the average of the two middle values
        train_median = (middle_1 + middle_2) / 2
        
        # Update the class with the largest of the two middle values to the average
        larger_mid_class = stage1_df[stage1_df[1] == max(middle_1, middle_2)][0].iloc[0]
        adjusted_df.loc[adjusted_df[0] == larger_mid_class, 1] = train_median

    else:
        train_median = np.median(train_num)

    # Create a DataFrame containing median sample count
    median_anchor = adjusted_df[adjusted_df[1] == train_median][0]

    # Get the indices of samples with the median sample count in the training set
    median_index = (train_label == int(median_anchor.iloc[0])).nonzero(as_tuple=True)[0]

    if classify_dim % 2 == 0:
        sampled_indices = torch.randperm(len(median_index))[:int(train_median)]
        median_index = median_index[sampled_indices]

    # Separate major (more samples than median) and minor classes (less samples than median)
    train_major = adjusted_df[adjusted_df[1] > train_median]
    train_minor = adjusted_df[adjusted_df[1] < train_median]

    # Calculate the fold difference between the median and the minor classes
    anchor_fold = np.array((train_median) / (train_minor[:][1]))
    
    # Get the class indices of minor and major classes
    minor_cts = train_minor[0].to_numpy()
    major_cts = train_major[0].to_numpy()

    # Get the data and labels corresponding to the indices of samples with the median sample count
    anchor_data = train_data[median_index.tolist(), :]
    anchor_label = train_label[median_index.tolist()]
    
    # Initialize new_data and new_label with the anchor_data and anchor_label
    new_data = anchor_data.to(device)
    new_label = anchor_label.to(device)

    ### Randomly downsample the major cell types
    j = 0
    for major_ct in major_cts:
        # Get the number of samples for the current major cell type
        Sample_num = np.array(train_major[1])[j]

        # Create a list of indices from 0 to anchor_num - 1
        N = range(int(Sample_num))

        # Get the indices of the samples belonging to the current major cell type
        index = (train_label == major_ct).nonzero(as_tuple = True)[0]

        # Get the data and labels of the samples belonging to the current major cell type
        Major_data = train_data[index.tolist(), :]
        Major_label = train_label[index.tolist()]

        # Randomly select 'train_median' samples from the current major cell type
        ds_index = random.sample(N, len(median_index))

        # Downsample the data and labels using the random indices
        Major_data = Major_data[ds_index, :]
        Major_label = Major_label[ds_index]

        # Concatenate the downsampled data and labels to the new data and label tensors
        new_data = torch.cat((new_data, Major_data.to(device)), 0)
        new_label = torch.cat((new_label, Major_label.to(device)), 0)

        j = j + 1

        # new_label_names.extend([label_to_name_mapping[major_ct][0]] * len(Major_data))


    ### Augment the minor cell types
    j = 0
    for minor_ct in minor_cts:
        # Calculate the augmentation fold and the remaining number of cells to reach train_median
        aug_fold = int((anchor_fold[j]))

        remaining_cell = train_median - aug_fold * np.array(train_minor[1])[j]

        # Get the indices of the samples belonging to the current minor cell type
        index = (train_label == minor_ct).nonzero(as_tuple = True)[0]

        # Get the data and labels of the samples belonging to the current minor cell type
        Minor_data = train_data[index.tolist(), :]
        Minor_label = train_label[index.tolist()]

        # Create a custom dataset and dataloader for the current minor cell type
        Minor_transform_dataset = MyDataset(Minor_data, Minor_label)
        Minor_dl = DataLoader(Minor_transform_dataset, batch_size = args.batch_size, shuffle = True, num_workers = 0, drop_last = False)

        # Generate new data samples using the VAE model
        reconstructed_data, reconstructed_label, real_data = get_vae_simulated_data_from_sampling(model, Minor_dl)

        # Clip the reconstructed_data values to the range of real_data values
        reconstructed_data[reconstructed_data > torch.max(real_data)] = torch.max(real_data)
        reconstructed_data[reconstructed_data < torch.min(real_data)] = torch.min(real_data)
        reconstructed_data[torch.isnan(reconstructed_data)] = torch.max(real_data)

        # Concatenate the generated data and labels to the new data and label tensors
        new_data = torch.cat((new_data, reconstructed_data.to(device)), 0)
        new_label = torch.cat((new_label, reconstructed_label.to(device)), 0)

        # new_label_names.extend([label_to_name_mapping[minor_ct][0]] * len(reconstructed_data))

        # Repeat the data generation process for the remaining augmentation folds
        for i in range(aug_fold - 1):
            # Generate new data samples using the VAE model
            reconstructed_data, reconstructed_label, real_data = get_vae_simulated_data_from_sampling(model, Minor_dl)

            # Clip the reconstructed_data values to the range of real_data values
            reconstructed_data[reconstructed_data > torch.max(real_data)] = torch.max(real_data)
            reconstructed_data[reconstructed_data < torch.min(real_data)] = torch.min(real_data)
            reconstructed_data[torch.isnan(reconstructed_data)] = torch.max(real_data)

            # Concatenate the generated data and labels to the new data and label tensors
            new_data = torch.cat((new_data, reconstructed_data.to(device)), 0)
            new_label = torch.cat((new_label, reconstructed_label.to(device)), 0)

            # new_label_names.extend([label_to_name_mapping[minor_ct][0]] * len(reconstructed_data))

        # Generate new data samples for the remaining cells
        reconstructed_data, reconstructed_label, real_data = get_vae_simulated_data_from_sampling(model, Minor_dl)

        # Clip the reconstructed_data values to the range of real_data values
        reconstructed_data[reconstructed_data > torch.max(real_data)] = torch.max(real_data)
        reconstructed_data[reconstructed_data < torch.min(real_data)] = torch.min(real_data)
        reconstructed_data[torch.isnan(reconstructed_data)] = torch.max(real_data)

        # Add the remaining cells
        # Create a list of indices from 0 to the number of samples in the current minor cell type - 1
        N = range(int(np.array(train_minor[1])[j]))

        # Randomly select the remaining cells to reach train_median
        ds_index = random.sample(N, int(remaining_cell))

        # Downsample the reconstructed_data and reconstructed_label using the random indices
        reconstructed_data = reconstructed_data[ds_index, :]
        reconstructed_label = reconstructed_label[ds_index]

        # Concatenate the downsampled data and labels to the new data and label tensors
        new_data = torch.cat((new_data, reconstructed_data.to(device)), 0)
        new_label = torch.cat((new_label, reconstructed_label.to(device)), 0)

        # Final check
        curr_cell_count = (new_label == minor_ct).sum().item()
        deficit = int(len(median_index) - curr_cell_count)

        if deficit > 0:
            # Generate new data samples using the VAE model
            reconstructed_data, reconstructed_label, real_data = get_vae_simulated_data_from_sampling(model, Minor_dl)

            # Clip the reconstructed_data values to the range of real_data values
            reconstructed_data[reconstructed_data > torch.max(real_data)] = torch.max(real_data)
            reconstructed_data[reconstructed_data < torch.min(real_data)] = torch.min(real_data)
            reconstructed_data[torch.isnan(reconstructed_data)] = torch.max(real_data)

             # Concatenate the generated data and labels to the new data and label tensors
            new_data = torch.cat((new_data, reconstructed_data.to(device)), 0)
            new_label = torch.cat((new_label, reconstructed_label.to(device)), 0)

        # new_label_names.extend([label_to_name_mapping[minor_ct][0]] * len(reconstructed_data))

        j = j + 1

    new_data, new_label = clip_counts(new_data, new_label, len(median_index))

    return new_data, new_label, new_label_names
