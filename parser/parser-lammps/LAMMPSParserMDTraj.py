import os, sys
import mdtraj as md
# from LAMMPSParserInput import readDumpFileName



# fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName()  # recover traj file name



def MDTrajParser(fNameTraj):

    ### LOADING TRAJECTORY AND TOPOLOGY
    mdTrajectory =  md.load(os.path.dirname(os.path.abspath(sys.argv[1])) + '/' + fNameTraj, top='top.pdb')
    mdTopology = md.load_topology('top.pdb')
    ###

    os.remove('top.pdb')

    MDTrajAtomPosition = mdTrajectory.xyz
    MDTrajAtomPosition = MDTrajAtomPosition.tolist()

    MDTrajSimulationCell = mdTrajectory.unitcell_lengths
    MDTrajSimulationCell = MDTrajSimulationCell.tolist()

    print MDTrajAtomPosition[0][0]
    print len(MDTrajAtomPosition[0][0])
    print len(MDTrajAtomPosition), len(MDTrajSimulationCell)

    return (MDTrajAtomPosition, MDTrajSimulationCell)