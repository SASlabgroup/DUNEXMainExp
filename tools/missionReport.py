# deploymentReport.py
'''
@edwinrainville

Description: This script goes through all data after it has been offloaded and makes a daily report of the data

The daily report will include:
1. List of microSWIFTs that were deployed that day
2. A map all drifter tracks - color coded by time 

'''

# Import Statements
import pylatex
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import pandas
import microSWIFTTools
import datetime

def make_report(doc, mission_num):
    '''
    @edwinrainville

    Builds mission report from data and the metaData spreadsheet.
    '''

    # Add document header
    header = pylatex.PageStyle("header")
    # Create left header
    with header.create(pylatex.Head("L")):
        header.append("Date : {}".format(datetime.date.today()))
    # Create center header
    with header.create(pylatex.Head("C")):
        header.append('DUNEX Project - University of Washington Group')
    # Create right header
    with header.create(pylatex.Head("R")):
        header.append(pylatex.simple_page_number())
    # Add Document Preamble and formatting
    doc.preamble.append(header)
    doc.change_document_style("header")

    # Add Heading
    with doc.create(pylatex.MiniPage(align='c')):
        doc.append(pylatex.LargeText(pylatex.utils.bold('microSWIFT Mission Number {}'.format(mission_num))))

    # Add Overview to document
    with doc.create(pylatex.Section('Overview')):
        # Time info     
        start_time = 'XXX'
        end_time = 'XXX'
        doc.append('Mission number {0} was from {1} to {2}. '.format(mission_num, start_time, end_time))

        # Forecast info
        Hs_forecast = 'XXX'
        Tp_forecasted = 'XXX'
        Dp_forecasted = 'XXX'
        doc.append('The forecasted conditions for the day were a significant wave height of {0} meters, \
                    a peak period of {1} seconds and a peak direction of {2} degrees relative to true north. \
                    '.format(Hs_forecast, Tp_forecasted, Dp_forecasted))

        # Array type and deployment
        deployment = 'XXX'
        array = 'XXX'
        doc.append('Due to these forecasted conditions, we deployment the microSWIFTs via {0} in a \
                     {1} array. '.format(deployment, array))

        # People and Positions
        deployer_lead = 'XXX'
        retriever_lead = 'XXX'
        aide = 'XXX'
        notetaker = 'XXX'
        doc.append('The lead deployer during this mission was {0}, the lead retriever was {1}, the aide was {2} and the \
                    notetaker was {3}. '.format(deployer_lead, retriever_lead, aide, notetaker))
        
        # microSWIFT final check
        microSWIFT_check = 'XXX'
        doc.append('The final microSWIFT check that all lights were on and lids were secured before deployment was done by {}. '.format(microSWIFT_check))

        # microSWIFTs deployed
        microSWIFTs_deployed = 'XXX'
        total_num_microSWFITs = 'XXX'
        doc.append('The microSWIFTs deployed during this mission were {0} for a total of {1} drifters. '.format(microSWIFTs_deployed, total_num_microSWFITs))

        # Additional Notes
        additional_notes = 'XXX'
        doc.append('The additional notes during the deployment from {0} were the following: '.format(notetaker))
        with doc.create(pylatex.Itemize()) as itemize:
            itemize.add_item(additional_notes)
     
        # Data offload Notes
        data_offload_check = 'XXX'
        data_offload_notes = 'XXX'
        doc.append('The data from this mission was offloaded by {0} and the additional data offload notes were: '.format(data_offload_check))
        with doc.create(pylatex.Itemize()) as itemize:
            itemize.add_item(data_offload_notes)

    # Add Mission Map to the document
    with doc.create(pylatex.Section('Drift Track Map')):
        map = microSWIFTTools.mission_map()
        doc.append(map)

    # Add Mission Accelerations to the document
    with doc.create(pylatex.Section('Acceleration Time Series')):
        accels = microSWIFTTools.mission_accels()
        doc.append(accels)

if __name__ == '__main__':

    # Define project directory
    project_dir = '../'

    # Define data directory 
    data_dir = 'microSWIFT_data/'

    # Define mission number 
    mission_num = int(input('Enter mission Number: '))

    # Define mission directroy
    mission_dir = '/mission_{}/'.format(mission_num)
    mission_dir_path = project_dir + data_dir + mission_dir

    # Create report and fill it with analysis of todays data
    doc_name_str = 'mission_{}_report'.format(mission_num)
    doc_path_str = mission_dir_path + doc_name_str
    geometry_options = {"margin": "0.7in"}
    doc = pylatex.Document(doc_path_str, geometry_options=geometry_options)
    make_report(doc, mission_num)

    # Generate the report
    doc.generate_tex()
    doc.generate_pdf(clean_tex=False)
