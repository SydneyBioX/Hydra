#!/usr/bin/python

##############################################

# Manoj M Wagle (USydney, MIT CSAIL)

##############################################


import torch, os, sys
import torch.nn as nn
from tqdm import tqdm
from torch.autograd import Variable
import torch.nn.functional as F
from collections import Counter
import numpy as np
from .util import AverageMeter,accuracy,save_checkpoint,CrossEntropyLabelSmooth,KL_loss

device = torch.device("cuda" if torch.cuda.is_available()
                      else "mps" if getattr(torch, 'has_mps', False)
                      else "cpu")


##############################################

class CurriculumSampler(torch.utils.data.Sampler):
    def __init__(self, labels, num_epochs):
        self.labels = labels
        self.num_epochs = num_epochs
        self.current_epoch = 0

    def __iter__(self):
        weights = self._compute_weights()
        total = len(self.labels)
        indices = np.random.choice(total, size=total, p=weights)
        return iter(indices)

    def _compute_weights(self):
        label_counts = Counter(self.labels)
        max_count = max(label_counts.values())
        curriculum_ratio = self.current_epoch / self.num_epochs
        weights = [1 / (label_counts[label] * (1 + curriculum_ratio)) for label in self.labels]

        # Normalize weights
        weight_sum = sum(weights)
        weights = [weight / weight_sum for weight in weights]
        return weights
    
    def update_epoch(self, epoch):
        self.current_epoch = epoch



##############################################

def train_model(model, train_dl, lr, epochs, classify_dim=17, save_path="", save_filename="", feature_num=10000, use_balancing=False):
    criterion = nn.MSELoss().cuda()
    criterion_smooth_cty = CrossEntropyLabelSmooth().cuda()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best_each_celltype_top1 = []
    best_each_celltype_num = []
    train_each_celltype_num = []
    
    for i in range(classify_dim):
        best_each_celltype_top1.append(0)
        best_each_celltype_num.append(0)
        train_each_celltype_num.append(0)

    best_loss = float('inf')
    no_improv = 0
    patience = 10

    # Curriculum learning sampler initialization
    if use_balancing:
        sampler = CurriculumSampler(train_dl.dataset.labels, epochs)
    else:
        sampler = None

    for epoch in tqdm(range(1, epochs + 1)):
        model.train()
        train_top1 = AverageMeter('Acc@1', ':6.2f')
        train_loss = 0.0
        train_batches = 0

        # Update the sampler for curriculum learning at the start of each epoch
        if use_balancing:
            sampler.update_epoch(epoch)

        for batch_sample in train_dl:
            optimizer.zero_grad()

            x = batch_sample['data'].to(device)
            x = Variable(x)
            x = torch.reshape(x, (x.size(0), -1))
            train_label = Variable(batch_sample['label']).to(device)

            encoder_out, x_cty, mu, var = model(x)
            
            # Calculating loss
            if isinstance(model, nn.DataParallel) and hasattr(model.module, 'classifiers'):
                x_cty = torch.mean(torch.stack([x_cty_i.to(device) for x_cty_i in x_cty]), dim=0)
            elif hasattr(model, 'classifiers'):
                x_cty = torch.mean(torch.stack([x_cty_i.to(device) for x_cty_i in x_cty]), dim=0)
            else:
                x_cty = x_cty.to(device)
                
            loss1 = criterion(encoder_out, x.to(device)) + 1 / feature_num * KL_loss(mu, var)
            loss2 = criterion_smooth_cty(x_cty, train_label.to(device))

            loss = 0.9 * loss1 + 0.1 * loss2
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            train_batches += 1
            
            train_pred1, = accuracy(x_cty, train_label, topk=(1,))
            train_top1.update(train_pred1[0], 1)

            if epoch == 1:
                for j in range(classify_dim):
                    if len(train_label[train_label == j]) != 0:
                        train_each_celltype_num[j] += len(train_label[train_label == j])

        average_train_loss = train_loss / train_batches
        
        # print(f"Epoch {epoch}/{epochs} - Train Loss: {average_train_loss:.4f}, Train Accuracy: {train_top1.avg:.4f}")
        
        if average_train_loss < best_loss:
            best_loss = average_train_loss
            no_improv = 0  # Reset the no improve counter
            model_to_save = model.module if isinstance(model, nn.DataParallel) else model
            torch.save({
                'epoch': epoch,
                'state_dict': model_to_save.state_dict(),
                'optimizer': optimizer.state_dict(),
            }, os.path.join(save_path, save_filename))
        else:
            no_improv += 1

        # If no improvement for 'patience' epochs, stop training
        if no_improv > patience:
            print("Early stopping at epoch:", epoch)
            break

    return model, best_loss, train_each_celltype_num
