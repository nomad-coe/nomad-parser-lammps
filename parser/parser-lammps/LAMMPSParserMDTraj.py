import os, sys
import mdtraj as md
from LAMMPSParserTraj import trajFileOpen

skipTraj = trajFileOpen()     # if False no trajectory info is available here

# the python module MDTraj is used to parse binary trajectory file formats (e.g., dcd, xtc)

def MDTrajParser(fNameTraj,trajDumpStyle):

    if skipTraj == False:

        skipMDtraj = False
        if fNameTraj and (trajDumpStyle == "atom" or trajDumpStyle == "xyz" or trajDumpStyle == "lammpstrj"):  ## these trajectory styles are handled by LAMMPSParserTraj.py
            skipMDtraj = True

        if fNameTraj and skipMDtraj == False:

            ### LOADING TRAJECTORY AND TOPOLOGY
            mdTrajectory =  md.load(os.path.dirname(os.path.abspath(sys.argv[1])) + '/' + fNameTraj, top='top.pdb') ### top.pdb is a pdb topology for the current system, created within LAMMPSParserData.py
            mdTopology = md.load_topology('top.pdb')
            ###

            MDTrajAtomPosition = mdTrajectory.xyz                  # atomic positions per frame
            MDTrajAtomPosition = MDTrajAtomPosition.tolist()

            MDTrajSimulationCell = mdTrajectory.unitcell_vectors   # simulation cell per frame
            MDTrajSimulationCell = MDTrajSimulationCell.tolist()

        else:
            MDTrajAtomPosition   = None
            MDTrajSimulationCell = None

        os.remove('top.pdb')  # remove the temporary .pdb topology (required to load a traj file through MDTraj)

    else:
        MDTrajSimulationCell = []
        MDTrajAtomPosition   = []
        os.remove('top.pdb')  # remove the temporary .pdb topology (required to load a traj file through MDTraj)

    return (MDTrajAtomPosition, MDTrajSimulationCell)