#!/usr/bin/env python


### Import libraries
import pandas as pd
import glob
import numpy as np
import plotly.express as px
from matplotlib import pyplot as plt

### Define functions
def bin_image(image_matrix, bin_size):
    h, w = image_matrix.shape
    
    # Ensure dimensions are divisible by bin_size
    new_h = (h // bin_size) * bin_size
    new_w = (w // bin_size) * bin_size
    image_matrix = image_matrix[:new_h, :new_w]  # Crop excess pixels

    # Reshape and bin
    binned_matrix = image_matrix.reshape(new_h // bin_size, bin_size, new_w // bin_size, bin_size).mean(axis=(1, 3))
    return binned_matrix

def make_matrix(column):
    image = np.zeros(shape=(X,Y))

    for element in range(0, len(df[column])):
        for x in range(df['XMin'][element]-df['XMin'].min(), df['XMax'][element]-df['XMax'].min()):
            for y in range(df['YMin'][element]-df['YMin'].min(), df['YMax'][element]-df['YMax'].min()):
                image[x, y] = df[column][element]*255
    return image

# Fetch directories
paths = glob.glob("../../raw_data/HALO_outputs/*")

# Make figures for masking and layer partitioning 
masking_metric_list = ['DAPI (DAPI) Nucleus Intensity','MBP (Opal 520) Cytoplasm Intensity']

for n in range(0,len(paths)):
    for masking_metric in masking_metric_list: 
        bin_size = 128
        df = pd.read_csv(paths[n])
        sample_id = paths[n].split('Br')[1].split('_DLPFC')[0]

        # Determine number of elements in X of matrix
        X = df['XMax'].max()-df['XMin'].min()

        # Determin number of elements in Y of matrix
        Y = df['YMax'].max()-df['YMin'].min()

        def make_matrix(column):
            image = np.zeros(shape=(X,Y))

            for element in range(0, len(df[column])):
                for x in range(df['XMin'][element]-df['XMin'].min(), df['XMax'][element]-df['XMax'].min()):
                    for y in range(df['YMin'][element]-df['YMin'].min(), df['YMax'][element]-df['YMax'].min()):
                        image[x, y] = df[column][element]*255
            return image

        image = make_matrix(masking_metric)

        pd.DataFrame(bin_image(image, bin_size = bin_size)).to_csv('../../processed_data/08_validation_oligos/Br%s_array_%s_%sx%sbinned_for_masking.txt' % (sample_id, masking_metric, bin_size, bin_size))

        from matplotlib import pyplot as plt

        plt.figure(figsize=(15, 15), dpi=300)
        plt.imshow(bin_image(image, bin_size = bin_size), interpolation = 'nearest', cmap = 'Greys')
        plt.gca().set_axis_off()
        plt.savefig('../../plots/08_validation_oligos/Br%s_array_DAPI_%sx%sbinned_for_masking.png' % (sample_id, bin_size, bin_size), bbox_inches = 'tight', dpi = 300)
        plt.show()
    
# make figures for actual data collection
metric_list = ['MBP_Oligo', # oligo counts
               'MBP (Opal 520) Nucleus Intensity', # MBP quantity
               'MBP (Opal 520) Cytoplasm Intensity', # MBP quantity
               'HPCAL1_L2_neuron', # Is HPCAL1 neuron
               'DAPI (DAPI) Positive'] # Is cell

for n in range(0,len(paths)):
    for m in range(0,len(metric_list)):
        
        bin_size = 1
        df = pd.read_csv(paths[n])
        sample_id = paths[n].split('Br')[1].split('_DLPFC')[0]
        metric = metric_list[m]

        # Determine number of elements in X of matrix
        X = df['XMax'].max()-df['XMin'].min()

        # Determin number of elements in Y of matrix
        Y = df['YMax'].max()-df['YMin'].min()

        image = make_matrix(metric)

        pd.DataFrame(bin_image(image, bin_size = bin_size)).to_csv('../../processed_data/08_validation_oligos/Br%s_array_%sx%sbinned_%s.txt' % (sample_id, bin_size, bin_size, metric))

        plt.figure(figsize=(15, 15), dpi=300)
        plt.imshow(bin_image(image, bin_size = bin_size), interpolation = 'nearest', cmap = 'Greys')
        plt.gca().set_axis_off()
        plt.savefig('../../plots/08_validation_oligos/Br%s_array_%sx%sbinned_%s.png' % (sample_id, bin_size, bin_size, metric), bbox_inches = 'tight', dpi = 300)
        plt.show()
        
        df[metric].hist(bins=20)
        plt.savefig('../../plots/08_validation_oligos/Br%s_histogram_%s.png' % (sample_id, metric), bbox_inches = 'tight', dpi = 300)
        plt.show()
        
