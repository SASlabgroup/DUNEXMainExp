# Import statements
import subprocess

# pullDUNEXData.py 
def main():
    '''
    @edwinrainville

    Description: Pulls microSWIFT data from the backed up copy to your local machine and updates anything on your local machine

    '''
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
        want_to_backup = input('Confirm you want to sync this data to your local machine(yes/no): ')

        if want_to_backup == 'yes':
            print('Data sync is starting')
            answer_val = 0

        elif want_to_backup == 'no':
            print('Exiting data sync')
            return

        else: 
            print('Enter either yes or no')

    # Define path to data directory on local machine
    microSWIFT_data_dir = '../microSWIFT_data/'

    # Define backup directory on bigwaves RAID
    bigwaves = 'bigwaves.apl.washington.edu'
    backup_dir = '/Volumes/Data/DuckFRF/DUNEXMainExp_2021/microSWIFT_data/'

    # Starting actual data backup
    print('Syncing data in {0} on bigwaves to {1} on local machine'.format(backup_dir, microSWIFT_data_dir))

    # rsync directory from remote backup to local machine
    # To download on mac OS use the command: brew install hudochenkov/sshpass/sshpass
    backup_process = subprocess.run(['sshpass', '-p', password, 'rsync', '-avz','{0}@{1}:{2}'.format(username, bigwaves, backup_dir), microSWIFT_data_dir]) 
    backup_process_rc = backup_process.returncode
    if backup_process_rc == 0:
        print('Data was successfully synced from bigwaves to local machine')
    else:
        print('Data was not successfully synced')

if __name__ == '__main__':
    main()