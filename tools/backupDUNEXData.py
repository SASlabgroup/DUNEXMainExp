# Import Statements
import subprocess
import logging
import datetime

def main():
    '''
    @edwinrainville

    Description: This script backs-up all DUNEXdata from the local machine to the remote bigwaves machine so that the data is 
    safe and secure.
    '''
    # Set up logging
    log_name = '../microSwIFT_data/data_backup.log'
    logging.basicConfig(filename=log_name, encoding='utf-8', level=logging.DEBUG)

    # Enter Login credentials for bigwaves RAID
    username = input('Enter bigwaves username: ')
    password = input('Enter bigwaves password: ')

    # Confirm that you are on the APL VPN
    answer_val = 1
    while answer_val != 0:
        want_to_backup = input('Confirm you are on the APL VPN(yes/no): ')

        if want_to_backup == 'yes':
            answer_val = 0

        elif want_to_backup == 'no':
            print('Exiting data backup')
            return

        else: 
            print('Enter either yes or no')

    # Confirm that you want to back up data
    answer_val = 1
    while answer_val != 0:
        want_to_backup = input('Confirm you want to back up this data(yes/no): ')

        if want_to_backup == 'yes':
            print('Data backup is starting')
            answer_val = 0

        elif want_to_backup == 'no':
            print('Exiting data backup')
            return

        else: 
            print('Enter either yes or no')

    # Define path to data directory
    microSWIFT_data_dir = '../microSWIFT_data/'

    # Define backup directory on bigwaves RAID
    bigwaves = 'bigwaves.apl.washington.edu'
    backup_dir = '/Volumes/Data/DuckFRF/DUNEXMainExp_2021/microSWIFT_data/'

    # Starting actual data backup
    logging.info('-------------- Data Backup at {} --------------'.format(datetime.datetime.utcnow()))
    logging.info('Data is being backed up by {}'.format(username))
    print('Backing up data in {0} on local machine to {1} on bigwaves'.format(microSWIFT_data_dir, backup_dir))

    # rsync directory
    backup_process = subprocess.run(['sshpass', '-p', password, 'rsync', '-avz', '--exclude', '*.nc', '--log-file={}'.format(log_name), microSWIFT_data_dir,'{0}@{1}:{2}'.format(username, bigwaves, backup_dir)]) 
    backup_process_rc = backup_process.returncode
    if backup_process_rc == 0:
        print('Data was successfully synced from local machine to bigwaves')
        logging.info('Data was successfully synced from local machine to bigwaves')
    else:
        print('Data was not successfully synced')
        logging.info('Data was not successfully synced')

# Run this main function as a script
if __name__ == '__main__':
    main()