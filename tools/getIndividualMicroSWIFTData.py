def main():
    '''
    @edwinrainville

    Description: This script will offload an individual microSWIFT from a mission and is used incase one microSWIFT 
    doesn't offload correctly when the main mission offloading occurs.
    '''
    # Import Statements
    import subprocess # Subprocess example: subprocess.run(["ls", "-l"])
    import pandas as pd
    import datetime
    import logging
    import datetime

    # Function to convert microSWIFT file name to datetime
    def get_microSWIFT_file_time(fname):
        import datetime

        # Convert Month string to month num
        month_str = fname[-21:-18]
        
        # January
        if month_str == 'Jan':
            month_num = '01'
        # February
        if month_str == 'Feb':
            month_num = '02'
        # March
        if month_str == 'Mar':
            month_num = '03'
        # April
        if month_str == 'Apr':
            month_num = '04'
        # May
        if month_str == 'May':
            month_num = '05'
        # June
        if month_str == 'Jun':
            month_num = '06'
        # July
        if month_str == 'Jul':
            month_num = '07'
        # August
        if month_str == 'Aug':
            month_num = '08'
        # September
        if month_str == 'Sep':
            month_num = '09'
        # October 
        if month_str == 'Oct':
            month_num = '10'
        # November
        if month_str == 'Nov':
            month_num = '11'
        # December
        if month_str == 'Dec':
            month_num = '12'

        # Compute Datetime
        date_str = '{0}-{1}-{2}T{3}:{4}:{5}'.format(fname[-18:-14], month_num, fname[-23:-21], fname[-13:-11], fname[-11:-9], fname[-9:-7])
        microSWIFT_file_time = datetime.datetime.fromisoformat(date_str)
        return microSWIFT_file_time

    # Define microSWIFT IP address and Password
    IP="192.168.0."
    PASSWORD="1013ne40th"

    # Define record Window Length 
    record_window_length = datetime.timedelta(hours=1)

    # Define project directory
    project_dir = '../'

    # Define Metadata Excel sheet name
    metadata_name = 'DUNEXMainExp_notes.xlsx'

    # Combine file name and project Directory
    metadata_filename = project_dir + metadata_name

    # User input for mission number 
    mission_num = int(input('Enter Mission Number: '))

    # User input for microSWIFT number to offload
    microSWIFT_num = int(input('Enter microSWIFT Number: '))

    # Create Data Directory for Mission
    mission_dir_str =  project_dir + "/microSWIFT_data/mission_{}".format(mission_num)
    subprocess.run(["mkdir", "-p", mission_dir_str])

    # Set up Data Offload Logging file
    log_name = mission_dir_str + '/data_offload.log'
    logging.basicConfig(filename=log_name, encoding='utf-8', level=logging.DEBUG)
    logging.info('------------ Mission {0} - microSWIFT {1} Individual Data Offload ------------'.format(mission_num, microSWIFT_num))

    # Create dataframe object from DUNEX MetaData SpreadSheet
    dunex_xlsx = pd.read_excel(metadata_filename)

    # Read in Start Time and convert to datetime
    start_time = datetime.datetime.fromisoformat(dunex_xlsx['Start Time'].iloc[mission_num])
    logging.info('Mission Start Time is: {}'.format(start_time))

    # Read in End Time and convert to datetime
    end_time = datetime.datetime.fromisoformat(dunex_xlsx['End Time'].iloc[mission_num])
    logging.info('Mission End Time is: {}'.format(end_time))

    # Ping microSWIFT to see if it is on the network
    microSWIFT_ip_address = IP + str(microSWIFT_num)
    ping = subprocess.run(['ping', '-c', '2', microSWIFT_ip_address])
    ping_val = ping.returncode

    # If microSWIFT is on network (return code from process is zero)
    if ping_val == 0:
        logging.info('microSWIFT {} is online'.format(microSWIFT_num))

        # Make Directory for this microSWIFT
        microSWIFT_dir_str =  mission_dir_str + '/microSWIFT_{}'.format(microSWIFT_num)
        subprocess.run(["mkdir", "-p", microSWIFT_dir_str])

        # Copy microSWIFT log into the microSWIFT directory
        # To download on mac OS use the command: brew install hudochenkov/sshpass/sshpass
        log_offload_process = subprocess.run(['sshpass', '-p', PASSWORD, 'scp', 'pi@{}:/home/pi/microSWIFT/logs/microSWIFT.log'.format(microSWIFT_ip_address), microSWIFT_dir_str ]) 
        log_offload_process_rc = log_offload_process.returncode
        if log_offload_process_rc == 0:
            logging.info('--- microSWIFT.log offloaded')
        else:
            logging.info('--- microSWIFT.log could not be offloaded')

        # Get list of all data files on microSWIFT
        list_of_data_files = subprocess.run(['sshpass', '-p', PASSWORD, 'ssh', 'pi@{}'.format(microSWIFT_ip_address), 'ls ~/microSWIFT/data'], stdout=subprocess.PIPE, text=True).stdout.splitlines()

        # Sort through each file to see if it is within the mission (within 1 record window)
        for file_name in list_of_data_files:
            # Get aware time object of when file was created
            file_time = get_microSWIFT_file_time(file_name)

            # Compare to see if it was within the mission time frame or within one burst length of the mission start time 
            if (file_time >= (start_time - record_window_length)  and file_time <= end_time):
                subprocess.run(['sshpass', '-p', PASSWORD, 'scp', 'pi@{0}:/home/pi/microSWIFT/data/{1}'.format(microSWIFT_ip_address, file_name), microSWIFT_dir_str])
                logging.info('--- {0} was copied to microSWIFT {1} data directory'.format(file_name, microSWIFT_num))

            else:
                continue

    else:
        logging.info('microSWIFT {} is offline'.format(microSWIFT_num))

if __name__ == "__main__":
    main()