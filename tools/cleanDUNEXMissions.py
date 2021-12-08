# Clean each DUNEX mission
import numpy as np
import logging
import netCDF4 as np

def main():
    '''
    @edwinrainville

    Description: Clean each DUNEX dataset by masking invalid points

    '''
    # Set up logger
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', encoding='utf-8', level=logging.DEBUG)

    # Define number of missions
    number_of_missions = 81

    # Load in each mission netCDF and apply masks 


if __name__ == '__main__':
    main()