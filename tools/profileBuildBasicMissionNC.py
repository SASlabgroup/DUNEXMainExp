# Profile buildBasicMissionNC.py
import cProfile,pstats
import sys

# Import DUNEX Tools
sys.path.append('..')
from tools import buildBasicMissionNC

def main():
    '''
    @edwinrainville
    
    '''
    profiler = cProfile.Profile()
    profiler.enable()
    mission_num = 20
    buildBasicMissionNC.main(mission_num=mission_num)
    profiler.disable()
    with open("profilingStatsAsText.txt", "w") as f:
        ps = pstats.Stats(profiler, stream=f).sort_stats('tottime')
        ps.print_stats()

if __name__ == '__main__':
    main()