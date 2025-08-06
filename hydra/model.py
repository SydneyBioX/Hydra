#!/usr/bin/python

##############################################

# Manoj M Wagle (USydney, MIT CSAIL)

##############################################


import torch
import torch.nn as nn
global mu
global var


##################################################

class Mish(nn.Module):
    def forward(self, x):
        return x * torch.tanh(torch.nn.functional.softplus(x))


##################################################

class LinBnDrop(nn.Sequential):
    """Module grouping `BatchNorm1d`, `Dropout` and `Linear` layers"""
    def __init__(self, n_in, n_out, bn=True, p=0., act=None, lin_first=True):
        self.p = p  # Store the dropout probability as an attribute
        layers = [nn.BatchNorm1d(n_out if lin_first else n_in)] if bn else []
        if p != 0: 
            layers.append(nn.Dropout(p))
        lin = [nn.Linear(n_in, n_out, bias=not bn)]
        if act is not None: lin.append(act)
        layers = lin+layers if lin_first else layers+lin
        super().__init__(*layers)


##################################################

class Encoder(nn.Module):
    """Encoder for 2 modalities (CITE-seq and SHARE-seq)"""
    def __init__(self, nfeatures_modality1=10703, nfeatures_modality2=192, hidden_modality1=185,  hidden_modality2=15, z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.nfeatures_modality2 = nfeatures_modality2
        self.encoder_modality1 = LinBnDrop(nfeatures_modality1, hidden_modality1, p=0.2, act=Mish())
        self.encoder_modality2 = LinBnDrop(nfeatures_modality2, hidden_modality2, p=0.2, act=Mish())
        self.encoder = LinBnDrop(hidden_modality1 + hidden_modality2, z_dim,  p=0.2, act=Mish())
        self.weights_modality1 = nn.Parameter(torch.rand((1,nfeatures_modality1)) * 0.001, requires_grad=True)
        self.weights_modality2 = nn.Parameter(torch.rand((1,nfeatures_modality2)) * 0.001, requires_grad=True)

        self.fc_mu = LinBnDrop(z_dim,z_dim, p=0.2)
        self.fc_var = LinBnDrop(z_dim,z_dim, p=0.2)

        self.mu = None
        self.var = None
        
    def reparameterize(self, mu, logvar):
        device = mu.device
        std = torch.exp(0.5 * logvar).to(device)
        eps = torch.randn_like(std).to(device)
        return eps * std + mu
    
    def forward(self, x):
        x_modality1 = self.encoder_modality1(x[:, :self.nfeatures_modality1]*self.weights_modality1)
        x_modality2 = self.encoder_modality2(x[:, self.nfeatures_modality1:]*self.weights_modality2)
        x = torch.cat([x_modality1, x_modality2], 1)
        x = self.encoder(x)
        self.mu = self.fc_mu(x)
        self.var = self.fc_var(x)
        x = self.reparameterize(self.mu, self.var)
        return x
    

class Decoder(nn.Module):
    """Decoder for 2 modalities (CITE-seq and SHARE-seq)"""
    def __init__(self, nfeatures_modality1=10703, nfeatures_modality2=192,  hidden_modality1=185,  hidden_modality2=15, z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.nfeatures_modality2 = nfeatures_modality2
        self.decoder1 = LinBnDrop(z_dim, nfeatures_modality1, act=Mish())
        self.decoder2 = LinBnDrop(z_dim, nfeatures_modality2,  act=Mish())
     
    def forward(self, x):
        x_rna = self.decoder1(x)
        x_adt = self.decoder2(x)
        x = torch.cat((x_rna,x_adt),1)
        return x


class Autoencoder_CITEseq_Step1(nn.Module):
    def __init__(self, nfeatures_rna=0, nfeatures_adt=0,  hidden_rna=185,  hidden_adt=15, z_dim=20,classify_dim=17, num_classifiers=1):
        """Autoencoder for CITE-seq Step-2"""
        super().__init__()
        self.encoder = Encoder(nfeatures_rna, nfeatures_adt, hidden_rna,  hidden_adt, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder(nfeatures_rna, nfeatures_adt, hidden_rna,  hidden_adt, z_dim)
    
    def forward(self, x):
        x = self.encoder(x)
        x_cty = [classifier(x) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty,mu,var


class Autoencoder_CITEseq_Step2(nn.Module):
    def __init__(self, nfeatures_rna=0, nfeatures_adt=0,  hidden_rna=185,  hidden_adt=15, z_dim=20,classify_dim=17):
        """Autoencoder for CITE-seq Step-1"""
        super().__init__()
        self.encoder = Encoder(nfeatures_rna, nfeatures_adt, hidden_rna,  hidden_adt, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder(nfeatures_rna, nfeatures_adt, hidden_rna,  hidden_adt, z_dim)
    
    def forward(self, x):
        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty,mu,var


##################################################

class Autoencoder_SHAREseq_Step1(nn.Module):
    def __init__(self, nfeatures_rna=0, nfeatures_atac=0,  hidden_rna=185,  hidden_atac=15, z_dim=20,classify_dim=17, num_classifiers=1):
        """Autoencoder for SHARE-seq Step-1"""
        super().__init__()
        self.encoder = Encoder(nfeatures_rna, nfeatures_atac, hidden_rna,  hidden_atac, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder(nfeatures_rna, nfeatures_atac, hidden_rna,  hidden_atac, z_dim)
    
    def forward(self, x):

        x = self.encoder(x)
        x_cty = [classifier(x) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty,mu,var
    

class Autoencoder_SHAREseq_Step2(nn.Module):
    def __init__(self, nfeatures_rna=0, nfeatures_atac=0,  hidden_rna=185,  hidden_atac=15, z_dim=20,classify_dim=17):
        """Autoencoder for SHARE-seq Step-2"""
        super().__init__()
        self.encoder = Encoder(nfeatures_rna, nfeatures_atac, hidden_rna,  hidden_atac, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder(nfeatures_rna, nfeatures_atac, hidden_rna,  hidden_atac, z_dim)
    
    def forward(self, x):

        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty,mu,var
    

##################################################

class Encoder_TEAseq(nn.Module):
    """Encoder for TEA-seq"""
    def __init__(self, nfeatures_rna=10703, nfeatures_adt=192, nfeatures_atac=192, hidden_rna=185, hidden_adt=30,  hidden_atac=185, z_dim=128):
        super().__init__()
        self.nfeatures_rna = nfeatures_rna
        self.nfeatures_adt = nfeatures_adt
        self.nfeatures_atac = nfeatures_atac

        self.encoder_rna = LinBnDrop(nfeatures_rna, hidden_rna, p=0.2, act=Mish())
        self.encoder_adt = LinBnDrop(nfeatures_adt, hidden_adt, p=0.2, act=Mish())
        self.encoder_atac = LinBnDrop(nfeatures_atac, hidden_atac, p=0.2, act=Mish())

        self.encoder = LinBnDrop(hidden_rna + hidden_adt +  hidden_atac, z_dim,  p=0.2, act=Mish())

        self.weights_rna = nn.Parameter(torch.rand((1,nfeatures_rna)) * 0.001, requires_grad=True)
        self.weights_adt = nn.Parameter(torch.rand((1,nfeatures_adt)) * 0.001, requires_grad=True)
        self.weights_atac = nn.Parameter(torch.rand((1,nfeatures_atac)) * 0.001, requires_grad=True)

        self.fc_mu = LinBnDrop(z_dim,z_dim)
        self.fc_var = LinBnDrop(z_dim,z_dim)
        self.mu = None
        self.var = None
        
    def reparameterize(self, mu, logvar):
        device = mu.device
        std = torch.exp(0.5 * logvar).to(device)
        eps = torch.randn_like(std).to(device)
        return eps * std + mu

    def forward(self, x):
        x_rna = self.encoder_rna(x[:, :self.nfeatures_rna]*self.weights_rna)
        x_adt = self.encoder_adt(x[:, self.nfeatures_rna:(self.nfeatures_rna+ self.nfeatures_adt)]*self.weights_adt)
        x_atac = self.encoder_atac(x[:, (self.nfeatures_rna+ self.nfeatures_adt):]*self.weights_atac)
        x = torch.cat([x_rna, x_adt, x_atac], 1)
        x = self.encoder(x)
        self.mu = self.fc_mu(x)
        self.var = self.fc_var(x)
        x = self.reparameterize(self.mu, self.var)
        return x
    

class Decoder_TEAseq(nn.Module):
    """Decoder for TEA-seq"""
    def __init__(self, nfeatures_rna=10703, nfeatures_adt=192, nfeatures_atac=10000, hidden_rna=185, hidden_adt=30, hidden_atac=185, z_dim=100):
        super().__init__()
        self.nfeatures_rna = nfeatures_rna
        self.nfeatures_adt = nfeatures_adt
        self.nfeatures_atac = nfeatures_atac
        self.decoder1 = LinBnDrop(z_dim, nfeatures_rna,  act=Mish())
        self.decoder2 = LinBnDrop(z_dim, nfeatures_adt,  act=Mish())
        self.decoder3 = LinBnDrop(z_dim, nfeatures_atac, act=Mish())

    def forward(self, x):
        x_rna = self.decoder1(x)
        x_adt = self.decoder2(x)
        x_atac = self.decoder3(x)
        x = torch.cat((x_rna,x_adt, x_atac), 1)
        return x
    

class Autoencoder_TEAseq_Step1(nn.Module):
    def __init__(self, nfeatures_rna=10000, nfeatures_adt=30, nfeatures_atac=10000, hidden_rna=185, hidden_adt=30,  hidden_atac=185, z_dim=100,classify_dim=17, num_classifiers=1):
        """ Autoencoder for TEA-seq Step-1"""
        super().__init__()
        self.encoder = Encoder_TEAseq(nfeatures_rna, nfeatures_adt, nfeatures_atac, hidden_rna,  hidden_adt, hidden_atac, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder_TEAseq(nfeatures_rna, nfeatures_adt, nfeatures_atac, hidden_rna,  hidden_adt, hidden_atac, z_dim)

    def forward(self, x):
        x = self.encoder(x)
        x_cty = [classifier(x) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

class Autoencoder_TEAseq_Step2(nn.Module):
    def __init__(self, nfeatures_rna=10000, nfeatures_adt=30, nfeatures_atac=10000, hidden_rna=185, hidden_adt=30, hidden_atac=185, z_dim=100, classify_dim=17):
        """ Autoencoder for TEA-seq Step-2"""
        super().__init__()
        self.encoder = Encoder_TEAseq(nfeatures_rna, nfeatures_adt, nfeatures_atac, hidden_rna, hidden_adt, hidden_atac, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder_TEAseq(nfeatures_rna, nfeatures_adt, nfeatures_atac, hidden_rna, hidden_adt, hidden_atac, z_dim)

    def forward(self, x):
        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var


##################################################
        
class Encoder_RNA_seq(nn.Module):
    """Encoder for RNA-seq"""
    def __init__(self, nfeatures_modality1=10703, hidden_modality1=185,  z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.encoder_modality1 = LinBnDrop(nfeatures_modality1, hidden_modality1, p=0.2, act=Mish())
        self.encoder = LinBnDrop(hidden_modality1, z_dim,  p=0.2, act=Mish())
        self.weights_modality1 = nn.Parameter(torch.rand((1,nfeatures_modality1)) * 0.001, requires_grad=True)

        self.fc_mu = LinBnDrop(z_dim,z_dim, p=0.2)
        self.fc_var = LinBnDrop(z_dim,z_dim, p=0.2)

        self.mu = None
        self.var = None
        
    def reparameterize(self, mu, logvar):
        device = mu.device
        std = torch.exp(0.5 * logvar).to(device)
        eps = torch.randn_like(std).to(device)
        return eps * std + mu
    
    def forward(self, x):
        x = self.encoder_modality1(x*self.weights_modality1)
        x = self.encoder(x)
        self.mu = self.fc_mu(x)
        self.var = self.fc_var(x).to(self.mu.device)
        x = self.reparameterize(self.mu, self.var)
        return x
    

class Decoder_RNA_seq(nn.Module):
    """Decoder for RNA-seq """
    def __init__(self, nfeatures_modality1=10703, hidden_modality1=185, z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.decoder1 = LinBnDrop(z_dim, nfeatures_modality1, act=Mish())

    def forward(self, x):
        x = self.decoder1(x)
        return x

    
class Autoencoder_RNAseq_Step1(nn.Module):
    def __init__(self, nfeatures_rna=0, hidden_rna=185, z_dim=20,classify_dim=17, num_classifiers=1):
        """ Autoencoder for RNA-seq Step-1"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_rna, hidden_rna, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder_RNA_seq(nfeatures_rna, hidden_rna, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        device = x.device
        x_cty = [classifier(x.to(device)) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

class Autoencoder_RNAseq_Step2(nn.Module):
    def __init__(self, nfeatures_rna=0, hidden_rna=185, z_dim=20,classify_dim=17):
        """ Autoencoder for RNA-seq Step-2"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_rna, hidden_rna, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder_RNA_seq(nfeatures_rna, hidden_rna, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

##################################################
        
class Encoder_ADT_seq(nn.Module):
    """Encoder for ADT-seq"""
    def __init__(self, nfeatures_modality1=192, hidden_modality1=30,  z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.encoder_modality1 = LinBnDrop(nfeatures_modality1, hidden_modality1, p=0.2, act=Mish())
        self.encoder = LinBnDrop(hidden_modality1, z_dim,  p=0.2, act=Mish())
        self.weights_modality1 = nn.Parameter(torch.rand((1,nfeatures_modality1)) * 0.001, requires_grad=True)

        self.fc_mu = LinBnDrop(z_dim,z_dim, p=0.2)
        self.fc_var = LinBnDrop(z_dim,z_dim, p=0.2)

        self.mu = None
        self.var = None
        
    def reparameterize(self, mu, logvar):
        device = mu.device
        std = torch.exp(0.5 * logvar).to(device)
        eps = torch.randn_like(std).to(device)
        return eps * std + mu
    
    def forward(self, x):
        x = self.encoder_modality1(x*self.weights_modality1)
        x = self.encoder(x)
        self.mu = self.fc_mu(x)
        self.var = self.fc_var(x).to(self.mu.device)
        x = self.reparameterize(self.mu, self.var)
        return x
    

class Decoder_ADT_seq(nn.Module):
    """Decoder for ADT-seq """
    def __init__(self, nfeatures_modality1=192, hidden_modality1=30, z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.decoder1 = LinBnDrop(z_dim, nfeatures_modality1, act=Mish())

    def forward(self, x):
        x = self.decoder1(x)
        return x

    
class Autoencoder_ADTseq_Step1(nn.Module):
    def __init__(self, nfeatures_adt=0, hidden_adt=30, z_dim=20, classify_dim=17, num_classifiers=1):
        """ Autoencoder for ADT-seq Step-1"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_adt, hidden_adt, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder_RNA_seq(nfeatures_adt, hidden_adt, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        device = x.device
        x_cty = [classifier(x.to(device)) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

class Autoencoder_ADTseq_Step2(nn.Module):
    def __init__(self, nfeatures_adt=0, hidden_adt=30, z_dim=20,classify_dim=17):
        """ Autoencoder for ADT-seq Step-2"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_adt, hidden_adt, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder_RNA_seq(nfeatures_adt, hidden_adt, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

##################################################
        
class Encoder_ATAC_seq(nn.Module):
    """Encoder for ATAC-seq"""
    def __init__(self, nfeatures_modality1=192, hidden_modality1=185,  z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.encoder_modality1 = LinBnDrop(nfeatures_modality1, hidden_modality1, p=0.2, act=Mish())
        self.encoder = LinBnDrop(hidden_modality1, z_dim,  p=0.2, act=Mish())
        self.weights_modality1 = nn.Parameter(torch.rand((1,nfeatures_modality1)) * 0.001, requires_grad=True)

        self.fc_mu = LinBnDrop(z_dim,z_dim, p=0.2)
        self.fc_var = LinBnDrop(z_dim,z_dim, p=0.2)

        self.mu = None
        self.var = None
        
    def reparameterize(self, mu, logvar):
        device = mu.device
        std = torch.exp(0.5 * logvar).to(device)
        eps = torch.randn_like(std).to(device)
        return eps * std + mu
    
    def forward(self, x):
        x = self.encoder_modality1(x*self.weights_modality1)
        x = self.encoder(x)
        self.mu = self.fc_mu(x)
        self.var = self.fc_var(x).to(self.mu.device)
        x = self.reparameterize(self.mu, self.var)
        return x
    

class Decoder_ATAC_seq(nn.Module):
    """Decoder for ATAC-seq """
    def __init__(self, nfeatures_modality1=192, hidden_modality1=185, z_dim=128):
        super().__init__()
        self.nfeatures_modality1 = nfeatures_modality1
        self.decoder1 = LinBnDrop(z_dim, nfeatures_modality1, act=Mish())

    def forward(self, x):
        x = self.decoder1(x)
        return x

    
class Autoencoder_ATACseq_Step1(nn.Module):
    def __init__(self, nfeatures_atac=0, hidden_atac=185, z_dim=20, classify_dim=17, num_classifiers=1):
        """ Autoencoder for ATAC-seq Step-1"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_atac, hidden_atac, z_dim)
        self.classifiers = nn.ModuleList([nn.Linear(z_dim, classify_dim) for _ in range(num_classifiers)])
        self.decoder = Decoder_RNA_seq(nfeatures_atac, hidden_atac, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        device = x.device
        x_cty = [classifier(x.to(device)) for classifier in self.classifiers]
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
    

class Autoencoder_ATACseq_Step2(nn.Module):
    def __init__(self, nfeatures_atac=0, hidden_atac=185, z_dim=20,classify_dim=17):
        """ Autoencoder for ATAC-seq Step-2"""
        super().__init__()
        self.encoder = Encoder_RNA_seq(nfeatures_atac, hidden_atac, z_dim)
        self.classify = nn.Linear(z_dim, classify_dim)
        self.decoder = Decoder_RNA_seq(nfeatures_atac, hidden_atac, z_dim)
        
    def forward(self, x):
        x = self.encoder(x)
        x_cty = self.classify(x)
        x = self.decoder(x)
        mu = self.encoder.mu
        var = self.encoder.var
        return x, x_cty, mu, var
